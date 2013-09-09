from __future__ import division
import collections
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
from scipy.sparse import dia_matrix
import warnings
np.random.seed(seed=1)

# Supporting functions
gaussian = lambda z, height, position, hwhm: height * np.exp(-np.log(2) * ((z - position)/hwhm)**2)
H = lambda z: 0.5 * (1 - np.sign(z))
TH = lambda x, sigma, mu: np.where( x>(mu-sigma), 1, 0) * np.where(x<(mu+sigma), 1, 0)

def check_index_within_bounds(i, min_i, max_i):
    """Checks that the index specified (can be number or an iterable) is within the given range."""
    success = np.all((i>=min_i)*(i<=max_i))
    if success:
        return True
        
    if isinstance(i, collections.Iterable):
        # The index is array-like
        print "Index is out of bounds.\ni=%s" % i[np.where(np.logical_not((i>=min_i)*(i<=max_i)))]
    else:
        # The index is an number
        print "Index is out of bounds.\ni=%s" % i
    return False


class Mesh(object):
    """A 1D cell centered mesh defined by faces for the finite volume method."""
    
    def __init__(self, faces):
        super(Mesh, self).__init__()
        
        # Check for duplicated points
        if len(faces) != len(set(faces)):
            raise ValueError("The faces array contains duplicated positions. No cell can have zero volume so please update with unique face positions.")
        self.faces = np.array(faces)
        self.cells = 0.5 * (self.faces[0:-1] + self.faces[1:])
        self.J = len(self.cells)
        self.cell_widths = (self.faces[1:] - self.faces[0:-1])
    
    def h(self, i):
        """Returns the width of the cell at the specified index."""
        return self.cell_widths[i]
    
    def hm(self, i):
        """Distance between centroids in the backwards direction."""
        if not check_index_within_bounds(i,1,self.J-1):
            raise ValueError("hm index runs out of bounds")
        return (self.cells[i] - self.cells[i-1])
        
    def hp(self, i):
        """Distance between centroids in the forward direction."""
        if not check_index_within_bounds(i,0,self.J-2):
            raise ValueError("hp index runs out of bounds")
        return (self.cells[i+1] - self.cells[i])

class CellVariable(np.ndarray):
    """Representation of a variable defined at the cell centers. Provides interpolation functions to calculate the value at cell faces."""
    
    # http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    def __new__(cls, input_array, mesh=None):
        
        # If `input_array` is actually just a constant
        # convert it to an array of len the number of cells.
        try:
            len(input_array)
        except:
            input_array = input_array*np.ones(len(mesh.cells))
            
        obj = np.asarray(input_array).view(cls)
        obj.mesh = mesh
        return obj
    
    def __array_finalize__(self, obj):
        if obj is None: return
        self.mesh = getattr(obj, 'mesh', None)
        self.__get_items__ = getattr(obj, '__get_items__', None)
        
    def m(self, i):
        """Linear interpolation of the cell value at the right hand face i.e. along the _m_inus direction."""
        return self.mesh.h(i)/(2*self.mesh.hm(i))*self[i-1] + self.mesh.h(i-1)/(2*self.mesh.hm(i))*self[i]
    
    def p(self, i):
        """Linear interpolation of the cell value at the right hand face i.e. along the _p_lus direction."""
        return self.mesh.h(i+1)/(2*self.mesh.hp(i))*self[i] + self.mesh.h(i)/(2*self.mesh.hp(i))*self[i+1]

class AdvectionDiffusionModel(object):
    
    """A model for the advection-diffusion equation"""
    def __init__(self, faces, a, d, reaction=None, discretisation="central"):
        
        print "moo~"
        super(AdvectionDiffusionModel, self).__init__()
        
        self.mesh = Mesh(faces)
        self.d0 = CellVariable(d, mesh=self.mesh)           # Unmodifed diffusion coefficient
        self.a = CellVariable(a, mesh=self.mesh)
        self.discretisation = discretisation
        self.set_reaction_term(reaction)
        
        # Default value
        self.has_transient = True
        
        # Artificially modify the diffusion coefficient to introduce adpative discretisation
        self.kappa = self.kappa_for_adpative_upwinding()
        self.d = CellVariable(self.modified_diffusion_coefficient_for_adaptive_upwinding(self.kappa), mesh=self.mesh)
        print "Using kappa:\n", "min =", np.min(self.kappa), "   max =", np.max(self.kappa)
    
    def update_velocity(self, new_velocity):
        """When the velocity is updated we need to recalculate 
        the adaptive upwind scheme."""
        self.a = CellVariable(new_velocity, mesh=self.mesh)
        self.kappa = self.kappa_for_adpative_upwinding()
        self.d = CellVariable(self.modified_diffusion_coefficient_for_adaptive_upwinding(self.kappa), mesh=self.mesh)
        
    def kappa_for_adpative_upwinding(self):
        
        mu = self.peclet_number()
        #print "peclet number", mu
        if self.discretisation == "exponential":
            # We need to compute the numerator and denominator separately because
            # mathmatically we expect inf/inf = 1 (in the limit mu->oo), however
            # computationally inf/inf = nan
            numerator = np.exp(mu) + 1
            denominator = np.exp(mu) - 1
            ratio = numerator/denominator               # Here catch the cases where
            ratio[np.where(np.isinf(numerator))] = 1    # we had inf/inf and replaced
            ratio[np.where(np.isinf(denominator))] = 1  # with the limiting value.
            kappa = ratio - 2/mu;                       # One last problem to fix, here
            kappa[np.where(mu==0.0)] = 0.0              # catch divide by zero values.
            
        elif self.discretisation == "upwind":
            kappa_neg = np.where(self.a<0,-1,0)
            kappa_pos = np.where(self.a>0,1,0)
            kappa = kappa_neg + kappa_pos
        elif self.discretisation == "central":
            kappa = np.zeros(self.mesh.J)
        else:
            print "Please set `discretisation` to one of the following: `upwind`, `central` or `exponential`."
        #print kappa
        assert (np.min(kappa) >= -1.0) and (np.max(kappa) <= 1.0), "kappa for adaptive upwinding does not have a sensible value for all cells."
        return kappa
    
    def modified_diffusion_coefficient_for_adaptive_upwinding(self, kappa):
        return self.d0 + 0.5 * self.a * self.mesh.cell_widths * kappa
        
    def peclet_number(self):
        
        mu = self.a * self.mesh.cell_widths / self.d0
        if self.discretisation == "central":
            if np.max(np.abs(mu)) >= 1.5 and np.max(np.abs(mu)) < 2.0:
                warnings.warn("\n\nThe Peclet number is %g, this is getting close to the limit of mod 2.")
            elif np.max(np.abs(mu)) > 2:
                warnings.warn("\n\nThe Peclet number (%g) has exceeded the maximum value of mod 2 for the central discretisation scheme. Suggest use change to upwind or an exponentially fitted discretisation." % (np.max(mu),) )
            
        return mu
    
    def CFL_condition(self, tau):
        return self.a * tau/ self.mesh.cell_widths
    
    def set_reaction_term(self, reaction_term):
        """Warning the reaction term influences the boundary conditions values 
        make sure that this is set _before_ setting the boundary conditions."""
        self.reaction = np.zeros(self.mesh.J) if reaction_term is None else reaction_term
        assert len(self.reaction) == self.mesh.J
        
    def set_has_transient(self, requires_transient_term):
        """A flag that changes elements in the boundary conditions vecs/mats."""
        self.has_transient = requires_transient_term
        
    def set_boundary_conditions(self, left_flux=None, right_flux=None, left_value=None, right_value=None ):
        """Make sure this function is used sensibly otherwise the matrix will be ill posed."""
        
        self.left_flux = left_flux
        self.right_flux = right_flux
        self.left_value = left_value
        self.right_value = right_value
        
    def _interior_matrix_elements(self, i):
        # Interior coefficients for matrix equation
        ra = lambda i, a, d, m: 1./m.h(i)*(a.m(i)*m.h(i)/(2*m.hm(i)) + d.m(i)/m.hm(i))
        rb = lambda i, a, d, m: 1./m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
        rc = lambda i, a, d, m: 1./m.h(i)*(-a.p(i)*m.h(i)/(2*m.hp(i)) + d.p(i)/m.hp(i))
        return ra(i, self.a, self.d, self.mesh), rb(i, self.a, self.d, self.mesh), rc(i,self.a, self.d, self.mesh)
        
    def _robin_boundary_condition_matrix_elements_left(self):
        # Left hand side Robin boundary coefficients for matrix equation
        b1 = lambda a, d, m: 1./m.h(0)*(-a.p(0)*m.h(1)/(2*m.hp(0)) - d.p(0)/m.hp(0) )
        c1 = lambda a, d, m: 1./m.h(0)*(-a.p(0)*m.h(0)/(2*m.hp(0)) + d.p(0)/m.hp(0) )
        
        # Index and element value
        locations = [(0,0), (0,1)]
        values = ( b1(self.a, self.d, self.mesh ), 
                   c1(self.a, self.d, self.mesh ) )
        return tuple([list(x) for x in zip(locations, values)])
        
    def _robin_boundary_condition_matrix_elements_right(self, matrix=None):
        # Right hand side Robin boundary coefficients for matrix equation
        aJ = lambda a, d, m: 1./m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-1)/(2*m.hm(m.J-1)) + d.m(m.J-1)/m.hm(m.J-1) )
        bJ = lambda a, d, m: 1./m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-2)/(2*m.hm(m.J-1)) - d.m(m.J-1)/m.hm(m.J-1) )
        J = self.mesh.J  # Index and element value
        
        # Index and element value
        locations = [(J-1,J-2), (J-1,J-1)]
        values = ( aJ(self.a, self.d, self.mesh ), 
                   bJ(self.a, self.d, self.mesh ) )
        return tuple([list(x) for x in zip(locations, values)])
        
    def _robin_boundary_condition_vector_elements_left(self):
        # Index and boundary condition vector elements for Robin conditions
        location = [0]
        value = [self.left_flux/self.mesh.h(0)]
        return tuple([list(x) for x in zip(location, value)])
        
    def _robin_boundary_condition_vector_elements_right(self):
        # Index and boundary condition vector elements for Robin conditions
        location = [self.mesh.J-1]
        value = [-self.right_flux/self.mesh.h(self.mesh.J-1)]
        return tuple([list(x) for x in zip(location, value)])
        
    def _dirichlet_boundary_condition_matrix_elements_left(self):
        # Left hand side Robin boundary coefficients for matrix equation
        rb = lambda i, a, d, m: 1./m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
        rc = lambda i, a, d, m: 1./m.h(i)*(-a.p(i)*m.h(i)/(2*m.hp(i)) + d.p(i)/m.hp(i))
        
        # Index and element value
        locations = [(0,0), (0,1)]
        # values = ( rb(0, self.a, self.d, self.mesh ), 
        #            rc(0, self.a, self.d, self.mesh ) )
        values = (1, 0 ) #BUG: these where set to (0, 1)
        return tuple([list(x) for x in zip(locations, values)])
        
    def _dirichlet_boundary_condition_matrix_elements_right(self):
        # Right hand side Robin boundary coefficients for matrix equation
        ra = lambda i, a, d, m: 1./m.h(i)*(a.m(i)*m.h(i)/(2*m.hm(i)) + d.m(i)/m.hm(i))
        rb = lambda i, a, d, m: 1./m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
        J = self.mesh.J  # Index and element value
        
        # Index and element value
        locations = [(J-1,J-2), (J-1,J-1)]
        # values = ( ra(self.J-1, self.a, self.d, self.mesh ), 
        #            rb(self.J-1, self.a, self.d, self.mesh ) )
        values = ( 0, 1 )
        return tuple([list(x) for x in zip(locations, values)])
        
    
    def _dirichlet_boundary_condition_vector_elements_left(self):
        # Index and boundary condition vector elements for Dirichlet conditions
        # NB these are always zero, unless BCs are time varying
        location = [0]
        value =  [0] if self.has_transient else [-self.reaction[0] - self.left_value]
        return tuple([list(x) for x in zip(location, value)])
        
    def _dirichlet_boundary_condition_vector_elements_right(self):
        # Index and boundary condition vector elements for Dirichlet conditions
        # NB these are always zero, unless BCs are time varying
        location = [self.mesh.J-1]
        value = [0] if self.has_transient else [-self.reaction[-1] - self.right_value]
        return tuple([list(x) for x in zip(location, value)])
    
    def alpha_matrix(self):
        """The alpha matrix is used to mask boundary conditions values for Dirichlet
        conditions. Otherwise for a fully Neumann (or Robin) system it is equal to 
        the identity matrix."""
        a1 = 1 if self.left_value is None else 0 if self.has_transient else 1
        aJ = 1 if self.right_value is None else 0 if self.has_transient else 1
        diagonals = np.ones(self.mesh.J)
        diagonals[0]  = a1
        diagonals[-1] = aJ
        offsets = (0,)
        return sparse.diags((diagonals,), offsets)
    
    def beta_vector(self, column_vector=False):
        """Returns the robin boundary condition vector."""
        b = np.zeros(self.mesh.J)
        
        if self.left_flux is not None:
            left_bc_elements = self._robin_boundary_condition_vector_elements_left()
        
        if self.right_flux is not None:
            right_bc_elements = self._robin_boundary_condition_vector_elements_right()
        
        if self.left_value is not None:
            left_bc_elements = self._dirichlet_boundary_condition_vector_elements_left()
        
        if self.right_value is not None:
            right_bc_elements = self._dirichlet_boundary_condition_vector_elements_right()
        
        bcs = left_bc_elements + right_bc_elements
        for inx, value in bcs:
            b[inx] = value
        
        if column_vector:
            return np.matrix(b).reshape((self.mesh.J,1))
        return b
        
    def coefficient_matrix(self):
        """Returns the coefficient matrix which appears on the left hand side."""
        J = self.mesh.J
        m = self.mesh
        a = self.a
        d = self.d
        
        padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
        zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
        one  = np.array([1])    #
        
        if self.left_flux is not None:
            left_bc_elements = self._robin_boundary_condition_matrix_elements_left()
        
        if self.right_flux is not None:
            right_bc_elements = self._robin_boundary_condition_matrix_elements_right()
        
        if self.left_value is not None:
            left_bc_elements = self._dirichlet_boundary_condition_matrix_elements_left()
        
        if self.right_value is not None:
            right_bc_elements = self._dirichlet_boundary_condition_matrix_elements_right()
        
        # Use the functions to layout the matrix Note that the boundary
        # condition elements are set to zero, they are filled in as
        # the next step.
        inx = np.array(range(1,J-1))
        ra, rb, rc = self._interior_matrix_elements(inx)
        #                                 c1
        upper = np.concatenate([padding, zero, rc ]) 
        
        #                          b1           bJ
        central = np.concatenate([zero, rb, zero  ]) 
        
        #                               aJ
        lower = np.concatenate([ra, zero , padding])
        
        A = sparse.spdiags([lower, central, upper], [-1,0,1], J, J).todok()
        
        # Apply boundary conditions elements
        bcs = left_bc_elements + right_bc_elements
        for inx, value in bcs:
            A[inx] = value
        return dia_matrix(A)

class PoissonEquation(AdvectionDiffusionModel):
    """Poisson equation."""
    def __init__(self, faces, permittivity, charge):
        
        # For the Poisson equation the dielectic constant is the diffusion
        # coefficient, there is no advection term.
        a = 0
        d = permittivity
        self.permittivity = permittivity
        self.charge = charge
        super(PoissonEquation, self).__init__(faces, a, d, 
            discretisation="central")
            

def test_poisson_equation():
    
    cells = 200
    faces = np.linspace(0,2,cells+1)
    permittivity = 1
    electron_valency = -1
    electron_concentration = np.ones(cells)
    charge = electron_concentration * electron_valency
    
    model = PoissonEquation(faces, permittivity, charge)
    left_value = 0
    right_flux = 0
    model.set_boundary_conditions(left_value=left_value, right_flux=right_flux)
    model.set_has_transient(False)
    model.set_reaction_term(charge)
    
    M = model.coefficient_matrix()
    alpha = model.alpha_matrix().todok()
    beta = model.beta_vector(column_vector=True)
    
    # Modify for equation with out transient term
    #alpha[0,0] = 1
    #beta[0] = -charge[0] - left_value
    
    # alpha[-1,-1] = 1
    # beta[-1] = -charge[-1] - right_value
    
    rho = np.matrix(charge).reshape((cells,1))
    #V = linalg.spsolve( alpha*(permittivity * M) + beta, -charge)
    
    V = linalg.spsolve( alpha*permittivity*M, -np.asarray(alpha*rho) - np.asarray(beta) )
    
    import pylab
    x = model.mesh.cells
    pylab.figure(1)
    pylab.plot(x, V, "--", label="Numerical approximation")
    pylab.plot(x, x**2/2 - 2*x, "-", label="Analytical solution")
    pylab.legend()
    pylab.show()
    
    
def test_advection_diffusion_equation():
    
    def geo_series(n, r, min_spacing=0.01):
        total = 0
        series = []
        for i in range(n):
            
            if i == 0:
                total = 1
            else:
                total = total - total*r
            series.append(total)
            
        series = np.array(series)
        norm = series / (np.max(series) - np.min(series))
        series = norm - np.min(norm)
        series =  np.abs(series - 1)
        series_diff = np.gradient(series)
        inx = np.where(series_diff > min_spacing)
        series_diff[inx] = min_spacing
        series_reconstruct = np.cumsum(series_diff)
        if np.min(series_reconstruct) != 0.0:
            series_reconstruct = np.array([0] + series_reconstruct.tolist())
        if np.max(series_reconstruct) != 1.0:
            series_reconstruct = np.array(series_reconstruct.tolist() + [1])
        return series_reconstruct
        
    #faces = geo_series(200, 0.15)
    #print faces.shape, faces
    
    #faces = np.concatenate((np.array([-0.5]), np.sort(np.random.uniform(-0.5, 1, 50)), np.array([1])))
    faces = np.linspace(0, 1, 200)
    #faces = np.concatenate([np.linspace(0, 0.99, 50), np.logspace(np.log10(0.991), np.log10(1.0), 100)])
    mesh =  Mesh(faces)
    
    a = CellVariable(1, mesh=mesh)    # Advection velocity
    d = CellVariable(1e-3, mesh=mesh) # Diffusion coefficient
    k = 0.0005                         # Time step 
    theta = 0.5
    left_value = 1.0
    #left_flux = 0.0
    right_flux = 0.0
    
    # Initial conditions
    w_init = 0.5*TH(mesh.cells, 0.4, 0)
    w_init = np.sin(np.pi*mesh.cells)**100
    w_init[0] = left_value
    #w_init[0] = left_flux
    
    # Source term
    #s[int(np.median(range(mesh.J)))] = 0.0
    
    
    model = AdvectionDiffusionModel(faces, a, d, discretisation="exponential")
    model.set_boundary_conditions(left_value=1., right_value=0.)
    #model.set_boundary_conditions(left_flux=left_flux, right_flux=left_flux)
    M = model.coefficient_matrix()
    alpha = model.alpha_matrix()
    beta  = model.beta_vector()
    I = sparse.identity(model.mesh.J)
    
    # Construct linear system from discretised matrices, A.x = d
    A = I - k*theta*alpha*M
    d = (I + k*(1-theta)*alpha*M)*w_init + beta
    
    print "Peclet number", np.min(model.peclet_number()), np.max(model.peclet_number())
    #print "CFL condition", np.min(model.CFL_condition()), np.max(model.CFL_condition())
    
    # matplotlib for movie export
    # see, http://matplotlib.org/examples/animation/moviewriter.html
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.animation as manimation
    print manimation.writers.__dict__
    FFMpegWriter = manimation.writers['ffmpeg']
    metadata = dict(title='Movie Test', artist='Matplotlib', comment='Movie support!')
    writer = FFMpegWriter(fps=15, metadata=metadata)
    
    fig = plt.figure()
    l0, = plt.plot([],[], 'r-', lw=1)
    l1, = plt.plot([],[], 'k-o', markersize=4)
    plt.xlim(np.min(faces), np.max(faces))
    plt.ylim(0,1.2)
    l1.set_data(mesh.cells,w_init)
    
    # # Analytical solution for Dirichlet boundary conditions
    analytical_x = np.concatenate([np.array([np.min(faces)]), mesh.cells, np.array([np.max(faces)])])
    analytical_solution = np.concatenate([np.array([model.left_value]), (np.exp(a/d) - np.exp(mesh.cells*a/d))/(np.exp(a/d)-1), np.array([model.right_value]) ])
    #analytical_solution2 = np.concatenate([np.array([model.left_value]), (np.exp(a/model.d) - np.exp(mesh.cells*a/model.d))/(np.exp(a/model.d)-1), np.array([model.right_value]) ])
    
    w = w_init
    with writer.saving(fig, "fvm_advection_diffusion_1.mp4", 300):
        
        data = [mesh.cells]
        for i in range(2000):
            #w = linalg.spsolve(A.tocsc(), M * w + s)
            d = (I + k*(1-theta)*alpha*M)*w + beta
            w = linalg.spsolve(A, d)
            
            # # Change velocity field sign if half way... will probably break!
            # if i > 300:
            #     model.update_velocity(-10*a)
            #     M = model.coefficient_matrix()
            #     alpha = model.alpha_matrix()
            #     beta  = model.beta_vector()
            #     d = (I + k*(1-theta)*alpha*M)*w + beta
                
                
            if  i == 0:
                l1.set_data(mesh.cells,w_init)
                writer.grab_frame()
            
            if i % 10 == 0 or i == 0:
                l1.set_data(mesh.cells,w)
                #l0.set_data(analytical_x, analytical_solution)
                area = np.sum(w * mesh.cell_widths)
                print "#%d; t=%g; area=%g:" % (i, i*k,area)
                writer.grab_frame()
                data.append(np.asarray(w))
    
    import pylab
    slices = data[1::2]
    x = data[0]
    for i, slice_f in enumerate(slices):
        fraction =( i+1) / len(slices)
        print fraction
        pylab.plot(x, slice_f, "-o", color=pylab.cm.Blues(1-fraction))
    
    from matplotlib import rc
    rc('font',**{'family':'serif','serif':['Palatino']})
    rc('text', usetex=True)
    pylab.title("Time stepping the linear advection-diffusion equation")
    pylab.xlabel("x")
    pylab.ylabel("w(x,t)")
    pylab.ylim(ymin=0.0, ymax=1.2)
    pylab.savefig("linear_advection_diffusion_IMEX.pdf")
    pylab.show()
    
    np.savetxt("data.txt", np.asarray(np.matrix(data).transpose()))

if __name__ == '__main__':
    #test_poisson_equation()
    test_advection_diffusion_equation()
    pass

    