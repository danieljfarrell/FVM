from __future__ import division
import collections
import numpy as np
from scipy import sparse
from scipy.sparse import linalg
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
    def __init__(self, faces, a, d, k, theta=0.5, discretisation="central"):
        super(AdvectionDiffusionModel, self).__init__()
        
        self.mesh = Mesh(faces)
        self.a = CellVariable(a, mesh=self.mesh)
        self.d = CellVariable(d, mesh=self.mesh)
        self.k = k
        self.theta = theta
        self.discretisation = discretisation
        
        # Check Peclet number
        import warnings
        mu = self.peclet_number()
        if np.max(np.abs(mu)) >= 1.5 and np.max(np.abs(mu)) < 2.0:
            warnings.warn("\n\nThe Peclet number is %g, this is getting close to the limit of mod 2.")
        elif np.max(np.abs(mu)) > 2:
            warnings.warn("\n\nThe Peclet number (%g) has exceeded the maximum value of mod 2 for the central discretisation scheme." % (np.max(mu),) )
        
        # Check CFL condition
        CFL = self.CFL_condition()
        if np.max(np.abs(CFL)) > 0.5 and np.max(np.abs(CFL)) < 1.0:
            warnings.warn("\n\nThe CFL condition value is %g, it is getting close to the upper limit." % (np.max(CFL),) )
        elif np.max(np.abs(CFL)) > 1:
            warnings.warn("\n\nThe CFL condition value is %g, and has gone above the upper limit." % (np.max(CFL),) )
            
        if discretisation == "exponential":
            self.kappa = (np.exp(mu) + 1)/(np.exp(mu) - 1) - 2/mu;
            self.kappa[np.where(mu==0.0)] = 0
            self.kappa[np.where(np.isposinf(mu))] = 1
            self.kappa[np.where(np.isneginf(mu))] = -1
        elif discretisation == "upwind":
            kappa_neg = np.where(self.a<0,-1,0)
            kappa_pos = np.where(self.a>0,1,0)
            self.kappa = kappa_neg + kappa_pos
        elif discretisation == "central":
            self.kappa = np.zeros(self.mesh.J)
        else:
            print "Please set `discretisation` to one of the following: `upwind`, `central` or `exponential`."
        
        # Artificially modify the diffusion coefficient to introduce adpative discretisation
        self.d = self.d + 0.5 * self.a * self.mesh.cell_widths * self.kappa
        print "Using kappa", np.min(self.kappa), np.max(self.kappa)
        print self.kappa
    
    def peclet_number(self):
        return self.a * self.mesh.cell_widths / self.d
    
    def CFL_condition(self):
        return self.a * self.k / self.mesh.cell_widths
        
    def set_boundary_conditions(self, left_flux=None, right_flux=None, left_value=None, right_value=None ):
        """Make sure this function is used sensibly otherwise the matrix will be ill posed."""
        
        self.left_flux = left_flux
        self.right_flux = right_flux
        self.left_value = left_value
        self.right_value = right_value
        
    def _interior_elements(self, i):
        
        # Interior coefficients for matrix equation
        ra = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i)/(2*m.hm(i)) + d.m(i)/m.hm(i))
        rb = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
        rc = lambda i, a, d, m, k: k/m.h(i)*(-a.p(i)*m.h(i)/(2*m.hp(i)) + d.p(i)/m.hp(i))
        return ra(i, self.a, self.d, self.mesh, self.k), rb(i, self.a, self.d, self.mesh, self.k), rc(i,self.a, self.d, self.mesh, self.k)
        
    def _robin_boundary_condition_elements_left(self, matrix=None):
        
        # Left hand side Robin boundary coefficients for matrix equation
        b1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(1)/(2*m.hp(0)) - d.p(0)/m.hp(0) )
        c1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(0)/(2*m.hp(0)) + d.p(0)/m.hp(0) )
        
        # Index and element value
        if matrix == "A":
            return [(0,0),1-self.theta*b1(self.a, self.d, self.mesh, self.k)], [(0,1),-self.theta*c1(self.a, self.d, self.mesh, self.k)]
        elif matrix == "M":
            return [(0,0),1+(1-self.theta)*b1(self.a, self.d, self.mesh, self.k)], [(0,1),(1-self.theta)*c1(self.a, self.d, self.mesh, self.k)]
        else:
            raise ValueError("Please choose a valid option for the `matrix` keyword, you can choose from `A` or `M`.")
    
    def _robin_boundary_condition_elements_right(self, matrix=None):
        
        # Right hand side Robin boundary coefficients for matrix equation
        aJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-1)/(2*m.hm(m.J-1)) + d.m(m.J-1)/m.hm(m.J-1) )
        bJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-2)/(2*m.hm(m.J-1)) - d.m(m.J-1)/m.hm(m.J-1) )
        J = self.mesh.J  # Index and element value
        
        # Index and element value
        if matrix == "A":
            locations = [(J-1,J-2), (J-1,J-1)]
            values = [-self.theta*aJ(self.a, self.d, self.mesh, self.k), 1-self.theta*bJ(self.a, self.d, self.mesh, self.k)]
            return tuple([list(x) for x in zip(locations, values)])
        elif matrix == "M":
            locations = [(J-1,J-2), (J-1,J-1)]
            values = [(1-self.theta)*aJ(self.a, self.d, self.mesh, self.k), 1+(1-self.theta)*bJ(self.a, self.d, self.mesh, self.k)]
            return tuple([list(x) for x in zip(locations, values)])
        else:
            raise ValueError("A valid option for `matrix` is `A` or `M`.")
        
    
    def _robin_boundary_condition_vector_elements_left(self):
        # Index and boundary condition vector elements for Robin conditions
        location = [0]
        value = [self.k*self.left_flux/self.mesh.h(0)]
        return tuple([list(x) for x in zip(location, value)])
        
    def _robin_boundary_condition_vector_elements_right(self):
        # Index and boundary condition vector elements for Robin conditions
        location = [self.mesh.J-1]
        value = [-self.k*self.right_flux/self.mesh.h(self.mesh.J-1)] # BUG FIXED HERE, changed to -1*self.right_flux
        return tuple([list(x) for x in zip(location, value)])
    
    def _dirichlet_boundary_condition_elements_left(self, matrix=None):
        # Left hand side Dirichlet boundary coefficients for matrix equation
        if matrix == "A":
            locations = [(0,0), (0,1), (1,0)]
            values = [1,0,0]
            return tuple([list(x) for x in zip(locations, values)])
        elif matrix == "M":
            locations = [(0,0), (0,1), (1,0)]
            values = [0,0,0]
            return tuple([list(x) for x in zip(locations, values)])
        else:
            raise ValueError("Please choose a valid option for the `matrix` keyword, you can choose from `A` or `M`.")
            
    def _dirichlet_boundary_condition_elements_right(self, matrix=None):
        # Right hand side Dirichlet boundary coefficients for matrix equation
        J = self.mesh.J
        if matrix == "A":
            locations = [(J-2,J-1), (J-1,J-2), (J-1,J-1)]
            values = [0,0,1]
            return tuple([list(x) for x in zip(locations, values)])
        elif matrix == "M":
            locations = [(J-2,J-1), (J-2,J-1), (J-1,J-1)]
            values = [0,0,0]
            return tuple([list(x) for x in zip(locations, values)])
        else:
            raise ValueError("Please choose a valid option for the `matrix` keyword, you can choose from `A` or `M`.")
    
    def _dirichlet_boundary_condition_vector_elements_left(self):
        # Vector elements for Dirichlet boundary condition
        ra, rb, rc = self._interior_elements(1)
        locations = [0,1]
        values = [self.left_value, ra*self.left_value]
        return tuple([list(x) for x in zip(locations, values)])
        
    def _dirichlet_boundary_condition_vector_elements_right(self):
        # Vector elements for Dirichlet boundary condition
        J = self.mesh.J
        ra, rb, rc = self._interior_elements(J-2)
        locations = [J-2, J-1]
        values = [rc*self.right_value, self.right_value]
        return tuple([list(x) for x in zip(locations, values)])
        
    def A_matrix(self):
        """Returns the coefficient matrix which appears on the left hand side."""
        J = self.mesh.J
        k = self.k
        m = self.mesh
        a = self.a
        d = self.d
        t = self.theta
        
        padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
        zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
        one  = np.array([1])    #
        
        if self.left_flux is not None:
            left_bc_elements = self._robin_boundary_condition_elements_left(matrix="A")
        
        if self.right_flux is not None:
            right_bc_elements = self._robin_boundary_condition_elements_right(matrix="A")
        
        if self.left_value is not None:
            left_bc_elements = self._dirichlet_boundary_condition_elements_left(matrix="A")
        
        if self.right_value is not None:
            right_bc_elements = self._dirichlet_boundary_condition_elements_right(matrix="A")
        
        # Use the functions to layout the matrix
        inx = np.array(range(1,J-1))
        ra, rb, rc = self._interior_elements(inx)
        #                                 c1
        upper = np.concatenate([padding, zero, -t*rc ]) 
        
        #                          b1           bJ
        central = np.concatenate([zero, 1-t*rb, zero  ]) 
        
        #                               aJ
        lower = np.concatenate([-t*ra, zero , padding])
        
        A = sparse.spdiags([lower, central, upper], [-1,0,1], J, J).todok()
        
        bcs = left_bc_elements + right_bc_elements
        for inx, value in bcs:
            A[inx] = value
        
        return A
        
    def M_matrix(self):
        """Returns the coefficient matrix which appears on the right hand side."""
        J = self.mesh.J
        k = self.k
        m = self.mesh
        a = self.a
        d = self.d
        t = self.theta
        
        padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
        zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
        one  = np.array([1])    #
        
        if self.left_flux is not None:
            left_bc_elements = self._robin_boundary_condition_elements_left(matrix="M")
        
        if self.right_flux is not None:
            right_bc_elements = self._robin_boundary_condition_elements_right(matrix="M")
        
        if self.left_value is not None:
            left_bc_elements = self._dirichlet_boundary_condition_elements_left(matrix="M")
        
        if self.right_value is not None:
            right_bc_elements = self._dirichlet_boundary_condition_elements_right(matrix="M")
        
        # Use the functions to layout the matrix
        inx = np.array(range(1,J-1))
        ra, rb, rc = self._interior_elements(inx)
        #                                 c1
        upper = np.concatenate([padding, zero, (1-t)*rc ])
        
        #                          b1               bJ
        central = np.concatenate([zero, 1+(1-t)*rb, zero  ]) 
        
        #                                   aJ
        lower = np.concatenate([(1-t)*ra, zero , padding])
        
        A = sparse.spdiags([lower, central, upper], [-1,0,1], J, J).todok()
        
        bcs = left_bc_elements + right_bc_elements
        for inx, value in bcs:
            A[inx] = value
        return A
    
    def b_vector(self):
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
        return b


class PoissonModel(object):
    """A finite volume solution of Poission equation"""
    def __init__(self, faces, epsilon, rho):
        super(PoissonModel, self).__init__()
        self.mesh = Mesh(faces)
        self.epsilon = CellVariable(epsilon, mesh=self.mesh)
        self.rho = CellVariable(rho, mesh=self.mesh)
    
    def _interior_elements(self, i):
        
        # Interior coefficients for matrix equation
        la = lambda i, e, r, m:  1/m.h(i)*e.m(i)/m.hm(i)
        lb = lambda i, e, r, m: -1/m.h(i)*(e.m(i)/m.hm(i) + e.p(i)/m.hp(i))
        lc = lambda i, e, r, m:  1/m.h(i)*e.p(i)/m.hp(i)
        return la(i, self.epsilon, self.rho, self.mesh), lb(i, self.epsilon, self.rho, self.mesh), lc(i, self.epsilon, self.rho, self.mesh)
    
    def A_matrix(self):
        
        padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
        zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
        one  = np.array([1])    #
        
        # Use the functions to layout the matrix
        J = self.mesh.J
        inx = np.array(range(1,J-1))
        la, lb, lc = self._interior_elements(inx)
        upper = np.concatenate([padding, zero, lc ]) 
        central = np.concatenate([zero, lb, zero  ]) 
        lower = np.concatenate([la, zero , padding])
        A = sparse.spdiags([lower, central, upper], [-1,0,1], J, J).todok()
        
        # import pylab
        # pylab.matshow(A.todense())
        # pylab.show()
        # 
        # Dirichlet left
        A[0,0] = 1
        A[0,1] = 0
        # Neumann right
        
        
        la = 1/self.mesh.h(J-1)*self.epsilon.m(J-1)/self.mesh.hm(J-1)
        A[-1,-1] = -la
        A[-1,-2] =  la
        
        # pylab.matshow(A.todense())
        # pylab.show()
        return -1*A
        
    def d_vector(self):
        
        J = self.mesh.J
        inx = np.array(range(0,J))
        d = self.rho[inx]
        
        omega_L = 0
        sigma_R = 0
        d[0] = -omega_L # left Dirichlet
        d[-1] = 1/self.mesh.h(J-1)*self.epsilon[J-1]*sigma_R + self.rho[J-1] # right Neumann
        
        # import pylab
        # pylab.matshow(np.matrix(d))
        # pylab.show()
        # print d
        return d
        
if __name__ == '__main__':
    
    #faces = np.linspace(0-0.005, 2+0.005, 200)
    # faces = np.linspace(0, 2, 200)
    # pm = PoissonModel(faces, 1, -1)
    # A = pm.A_matrix()
    # d = pm.d_vector()
    # w = linalg.spsolve(A.tocsc(), d)
    # print w
    # import pylab
    # pylab.plot(pm.mesh.cells, w, "-b")
    # x = faces
    # pylab.plot(x, x**2/2 - 2*x, "-r" )
    # pylab.show()
    # x = pm.mesh.cells
    # error = w - (x**2/2 - 2*x)
    # print error
    # #error = (error / (x**2/2 - 2*x) - 1) * 100
    # pylab.plot(x, error)
    # pylab.show()
    
    faces = np.concatenate((np.array([-0.5]), np.sort(np.random.uniform(-0.5, 1, 50)), np.array([1])))
    #faces = np.linspace(0, 1, 50)
    mesh =  Mesh(faces)
    
    a = CellVariable(1, mesh=mesh) # Advection velocity
    d = CellVariable(1e-3, mesh=mesh) # Diffusion coefficient
    k = 0.01                             # Time step 
    
    model = AdvectionDiffusionModel(faces, a, d, k, discretisation="exponential")
    model.set_boundary_conditions(left_value=1., right_value=0.)
    #model.set_boundary_conditions(left_flux=0, right_flux=0)
    #model.set_boundary_conditions(left_value=1, right_flux=0)
    A = model.A_matrix()
    M = model.M_matrix()
    b = model.b_vector()
    
    print "Peclet number", np.min(model.peclet_number()), np.max(model.peclet_number())
    print "CFL condition", np.min(model.CFL_condition()), np.max(model.CFL_condition())
    
    # Initial conditions
    w_init = 0.5*TH(mesh.cells, 0.4, 0)
    w_init = np.sin(np.pi*mesh.cells)**100
    
    # Source term
    b[int(np.median(range(mesh.J)))] = 0.0
    
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
    
        for i in range(201):
            w = linalg.spsolve(A.tocsc(), M * w + b)
        
            if  i == 0:
                l1.set_data(mesh.cells,w_init)
                writer.grab_frame()
            
            if i %  1 == 0 or i == 0:
                l1.set_data(mesh.cells,w)
                #l0.set_data(analytical_x, analytical_solution)
                area = np.sum(w * mesh.cell_widths)
                print "#%d; t=%g; area=%g:" % (i, i*k,area)
                writer.grab_frame()



    