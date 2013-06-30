from __future__ import division
import collections
import numpy as np
#import pylab
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

def covert_images_to_animated_gif(target_dir, filename="animation2.gif"):
    from visvis.vvmovie.images2gif import writeGif
    from PIL import Image
    import os
    
    file_names = sorted((os.path.join(target_dir, fn) for fn in os.listdir(target_dir) if fn.endswith('.png')))
    # Numerical sort from file names
    numbers = [float(os.path.basename(item).split("_")[0]) for item in file_names]
    print numbers
    sorted_indexes = np.argsort(numbers)
    file_names = np.array(file_names)
    sorted_file_names = file_names[sorted_indexes]
    print sorted_file_names
    #images = [Image.open(fn) for fn in sorted_file_names] # MacOS sets a limit for the number of files that can be open, this will break sometimes.
    images = list()
    for fn in sorted_file_names:
        file_obj = open(fn, "rb")
        img = Image.open(file_obj)
        img.load()
        images.append(img)
        file_obj.close()
    #images.extend(reversed(images)) #infinit loop will go backwards and forwards.
    if len(images)==0:
        print "The image directory is empty. Could not create the gif."
    else:
        print "writing gif"
        writeGif(filename, images, duration=0.2)
    




left_face = lambda i, f, h, hm: h[i]/(2*hm(i))*f[i-1] + h[i-1]/(2*hm(i))*f[i] # Interpolation from cell to face values
right_face = lambda i, f, h, hp: h[i+1]/(2*hp(i))*f[i] + h[i]/(2*hp(i))*f[i+1]
def advection_diffusion_A_matrix(a, d, mesh, k, ra, rb, rc, left_corner, right_corner, J):
    padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
    zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
    one  = np.array([1])    #
    inx = np.array(range(1,J-2))
    b1, c1, a2 = [np.array([item]) for item in left_corner]
    c1J, aJ, bJ = [np.array([item]) for item in right_corner]
    upper = np.concatenate([padding, c1, -theta*rc(inx,a, d, mesh, k), c1J])
    inx = np.array(range(1,J-1))
    central = np.concatenate([b1, 1-theta*rb(inx,a, d, mesh, k), bJ])
    inx = np.array(range(2,J-1))
    lower = np.concatenate([a2, -theta*ra(inx,a, d, mesh, k), aJ, padding])
    return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

def advection_diffusion_M_matrix(a, d, mesh, k, ra, rb, rc, left_corner, right_corner, J):
    padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
    zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
    one  = np.array([1])    #
    inx = np.array(range(1,J-2))
    b1, c1, a2 = [np.array([item]) for item in left_corner]
    c1J, aJ, bJ = [np.array([item]) for item in right_corner]
    upper = np.concatenate([padding, c1, (1-theta)*rc(inx,a,d,mesh,k), c1J])
    inx = np.array(range(1,J-1))
    central = np.concatenate([b1, 1+(1-theta)*rb(inx,a,d,mesh,k),bJ])
    inx = np.array(range(2,J-1))
    lower = np.concatenate([a2, (1-theta)*ra(inx,a,d,mesh,k), aJ, padding])
    return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)





class Mesh(object):
    """A 1D cell centered mesh defined by faces for the finite volume method."""
    
    def __init__(self, faces):
        super(Mesh, self).__init__()
        self.faces = np.array(faces)
        self.cells = 0.5 * (self.faces[0:-1] + self.faces[1:])
        self.J = len(self.cells)
        self.cell_widths = self.faces[1:] - self.faces[0:-1]
    
    def h(self, i):
        """Returns the width of the cell at the specified index."""
        return self.cell_widths[i]
    
    def hm(self, i):
        """Distance between centroids in the backwards direction."""
        if not check_index_within_bounds(i,1,self.J-1):
            raise ValueError("hm index runs out of bounds")
        return self.cells[i] - self.cells[i-1]
        
    def hp(self, i):
        """Distance between centroids in the forward direction."""
        if not check_index_within_bounds(i,0,self.J-2):
            raise ValueError("hp index runs out of bounds")
        return self.cells[i+1] - self.cells[i]

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

class Model(object):
    
    """A model for the advection-diffusion equation"""
    def __init__(self, faces, a, d, k, theta=0.5, discretisation="central"):
        super(Model, self).__init__()
        
        self.mesh = Mesh(faces)
        self.a = CellVariable(a, mesh=mesh)
        self.d = CellVariable(d, mesh=mesh)
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
            self.kappa = np.ones(mesh.J)
        elif discretisation == "central":
            self.kappa = np.zeros(mesh.J)
        else:
            print "Please set `discretisation` to one of the following: `upwind`, `central` or `exponential`."
        
        # Artificially modify the diffusion coefficient to introduce adpative discretisation
        self.d = self.d + 0.5 * self.a * self.kappa
    
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
        
    def _interior_functions(self):
        
        # Interior coefficients for matrix equation
        ra = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i)/(2*m.hm(i)) + d.m(i)/m.hm(i))
        rb = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
        rc = lambda i, a, d, m, k: k/m.h(i)*(-a.p(i)*m.h(i)/(2*m.hp(i)) + d.p(i)/m.hp(i))
        return ra, rb, rc
        
    def _robin_boundary_condition_elements_left(self):
        # Left hand side Robin boundary coefficients for matrix equation
        b1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(1)/(2*m.hp(0)) - d.p(0)/m.hp(0) )
        c1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(0)/(2*m.hp(0)) + d.p(0)/m.hp(0) )
        return b1(self.a, self.d, self.mesh, self.k), c1(self.a, self.d, self.mesh, self.k)
        
    def _robin_boundary_condition_elements_right(self):
        # Right hand side Robin boundary coefficients for matrix equation
        aJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-1)/(2*m.hm(m.J-1)) + d.m(m.J-1)/m.hm(m.J-1) )
        bJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-2)/(2*m.hm(m.J-1)) - d.m(m.J-1)/m.hm(m.J-1) )
        return aJ(self.a, self.d, self.mesh, self.k), bJ(self.a, self.d, self.mesh, self.k)
    
    def _robin_boundary_condition_vector_elements_left(self):
        # Vector elements for Robin boundary condition
        return self.k*self.left_flux/self.mesh.h(0), 0
        
    def _robin_boundary_condition_vector_elements_right(self):
        # Vector elements for Robin boundary condition
        return 0, self.k*self.right_flux/self.mesh.h(self.mesh.J-1)
    
    def _dirichlet_boundary_condition_elements_left(self):
        # Left hand side Dirichlet boundary coefficients for matrix equation
        b1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(1)/(2*m.hp(0)) - d.p(0)/m.hp(0) )
        c1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(0)/(2*m.hp(0)) + d.p(0)/m.hp(0) )
        return b1(self.a, self.d, self.mesh, self.k), c1(self.a, self.d, self.mesh, self.k)
        
    def _dirichlet_boundary_condition_elements_right(self):
        # Right hand side Dirichlet boundary coefficients for matrix equation
        aJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-1)/(2*m.hm(m.J-1)) + d.m(m.J-1)/m.hm(m.J-1) )
        bJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-2)/(2*m.hm(m.J-1)) - d.m(m.J-1)/m.hm(m.J-1) )
        return aJ(self.a, self.d, self.mesh, self.k), bJ(self.a, self.d, self.mesh, self.k)
    
    def _dirichlet_boundary_condition_vector_elements_left(self):
        # Vector elements for Dirichlet boundary condition
        return self.k*self.left_flux/self.mesh.h(0), 0
        
    def _dirichlet_boundary_condition_vector_elements_right(self):
        # Vector elements for Dirichlet boundary condition
        return 0, self.k*self.right_flux/self.mesh.h(self.mesh.J-1)
        
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
        
        ra, rb, rc = self._interior_functions()
        
        if self.left_flux is not None:
            b1, c1 = self._robin_boundary_condition_elements_left()
        
        if self.right_flux is not None:
            aJ, bJ = self._robin_boundary_condition_elements_right()
        
        if self.left_value is not None:
            raise NotImplementedError()
        
        if self.right_value is not None:
            raise NotImplementedError()
        
        # Use the functions to layout the matrix
        inx = np.array(range(1,J-1))
        upper = np.concatenate([padding, np.array([-t*c1]), -t*rc(inx,a, d, mesh, k) ])
        inx = np.array(range(1,J-1))
        central = np.concatenate([np.array([1-t*b1]) , 1-t*rb(inx,a, d, mesh, k), np.array([1-t*bJ]) ])
        inx = np.array(range(1,J-1))
        lower = np.concatenate([-t*ra(inx,a, d, mesh, k), np.array([-t*aJ]) , padding])
        return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)
        
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
        
        ra, rb, rc = self._interior_functions()
        
        if self.left_flux is not None:
            b1, c1 = self._robin_boundary_condition_elements_left()
        
        if self.right_flux is not None:
            aJ, bJ = self._robin_boundary_condition_elements_right()
        
        if self.left_value is not None:
            raise NotImplementedError()
        
        if self.right_value is not None:
            raise NotImplementedError()
        
        # Use the functions to layout the matrix
        inx = np.array(range(1,J-1))
        upper = np.concatenate([padding, np.array([(1-t)*c1]), (1-t)*rc(inx,a, d, mesh, k) ])
        inx = np.array(range(1,J-1))
        central = np.concatenate([np.array([1+(1-t)*b1]) , 1+(1-t)*rb(inx,a, d, mesh, k), np.array([1+(1-t)*bJ]) ])
        inx = np.array(range(1,J-1))
        lower = np.concatenate([(1-t)*ra(inx,a, d, mesh, k), np.array([(1-t)*aJ]) , padding])
        return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)
    
    def b_vector(self):
        """Returns the robin boundary condition vector."""
        b = np.zeros(self.mesh.J)
        
        
        if self.left_flux is not None:
            left_bc = self.left_flux
            b1, b2 = self._robin_boundary_condition_vector_elements_left()
        
        if self.right_flux is not None:
            right_bc = self.left_flux
            b1J, bJ = self._robin_boundary_condition_vector_elements_right()
        
        if self.left_value is not None:
            raise NotImplementedError()
        
        if self.right_value is not None:
            raise NotImplementedError()
            
        b[0] = b1
        b[1] = b2
        b[-2] = b1J
        b[-1] = bJ
        return b

#faces = np.concatenate((np.array([-0.5]), np.sort(np.random.uniform(-0.5, 1, 50)), np.array([1])))
faces = np.linspace(-0.5, 1, 200)
mesh =  Mesh(faces)

#a = CellVariable(0.01*TH(mesh.cells, 0.1, 0), mesh=mesh)  # Advection velocity
a = CellVariable(0.01, mesh=mesh)
d = CellVariable(1e-4, mesh=mesh) # Diffusion coefficient
k = 0.05                               # Time step 

model = Model(faces, a, d, k, discretisation="exponential")
model.set_boundary_conditions(left_flux=0., right_flux=0.)
A = model.A_matrix()
M = model.M_matrix()
b = model.b_vector()

print "Peclet number", np.min(model.peclet_number()), np.max(model.peclet_number())
print "CFL condition", np.min(model.CFL_condition()), np.max(model.CFL_condition())
# Initial conditions
w_init = 0.5*TH(mesh.cells, 0.4, 0)

# Source term
#b[int(np.median(range(mesh.J)))] = 0.01

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
l2, = plt.plot([],[], 'k|')
l3, = plt.plot([],[], 'k|')

plt.xlim(np.min(faces), np.max(faces))
plt.ylim(-0.2,1.1)
l1.set_data(mesh.cells,w_init)

w = w_init
with writer.saving(fig, "writer_test.mp4", 300):
    
    for i in range(200):
        #b[-1] = w[-1] * -0.1
        #b[int(np.median(range(mesh.J)))] = np.random.uniform(0.1,0)
        w = linalg.spsolve(A.tocsc(), M * w + b)
        
        if  i == 0:
            l1.set_data(mesh.cells,w_init)
            writer.grab_frame()
            
        if i %  1 == 0 or i == 0:
            #pylab.ylim((0,0.07))
            #pylab.bar(faces[:-1], w, width=(faces[1:]-faces[:-1]), edgecolor='white')
            #pylab.plot(cells, w, "k", lw=3)
            l1.set_data(mesh.cells,w)
            #l2.set_data(faces[0:-1], w)
            #l3.set_data(faces[1:], w)
            #l0.set_data(analytical_x, analytical_solution)
            #pylab.plot(cells, a, "--g", lw=2)
            #pylab.plot(analytical_x, analytical_solution, "r-", lw=2)
            area = np.sum(w * mesh.cell_widths)
            print "#%d; area:" % (i,), area
            #fig.suptitle("Area: %g" % (area))
            writer.grab_frame()
            #pylab.savefig(os.path.join(working_directory, "%d_solution.png" % i), dpi=72)

# Write the animated gif to the script directory
#import os
#pwd = os.path.join(os.path.dirname(os.path.realpath(__file__)))
#covert_images_to_animated_gif("/tmp/output", filename=os.path.join(pwd, "solution.gif") )



    