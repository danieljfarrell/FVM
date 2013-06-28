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

mesh = Mesh([1,2,3,4,5,6])
print mesh.h(4)
print mesh.hp(3)
print mesh.hm(4)

a = CellVariable(1, mesh=mesh)
print a
print a.m(1)
print a.p(1)


faces = np.concatenate((np.array([-0.5]), np.sort(np.random.uniform(-0.5, 1, 50)), np.array([1])))
mesh =  Mesh(faces)
J = len(mesh.cells)

a = CellVariable(0.001, mesh=mesh)  # Advection velocity
d = CellVariable(0.0001, mesh=mesh) # Diffusion coefficient
k = 1           # Time step 
theta = 0.5     # Implicit/Explicit
assert(theta>=0 or theta<=1, "theta is out of sensible range.")


# Check Peclet number and CFL condition
mu = a * mesh.cell_widths / d
CFL = a * k / mesh.cell_widths
print "Peclet number", np.min(mu), np.max(mu)
print "CFL condition", np.min(CFL), np.max(CFL)


kappa = (np.exp(mu) + 1)/(np.exp(mu) - 1) - 2/mu;
kappa[np.where(mu==0.0)] = 0
kappa[np.where(np.isposinf(mu))] = 1
kappa[np.where(np.isneginf(mu))] = -1
d += 0.5 * a * mesh.cell_widths * kappa # Artificially alter the diffusion coefficient to bring in adaptive upwinding
assert(np.any(kappa<-1) and np.any(kappa>1), "kappa is out side of sensible to range.")

# Initial conditions
#w_init = gaussian(cells, 1, 0, 0.01)
w_init = 0.2*H(mesh.cells)


# Interior coefficients for matrix equation
ra = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i)/(2*m.hm(i)) + d.m(i)/m.hm(i))
rb = lambda i, a, d, m, k: k/m.h(i)*(a.m(i)*m.h(i-1)/(2*m.hm(i)) - a.p(i)*m.h(i+1)/(2*m.hp(i)) - d.m(i)/m.hm(i) - d.p(i)/m.hp(i))
rc = lambda i, a, d, m, k: k/m.h(i)*(-a.p(i)*m.h(i)/(2*m.hp(i)) + d.p(i)/m.hp(i))


# Left hand side Robin boundary coefficients for matrix equation
robin_b1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(1)/(2*m.hp(0)) - d.p(0)/m.hp(0) )
robin_c1 = lambda a, d, m, k: k/m.h(0)*(-a.p(0)*m.h(0)/(2*m.hp(0)) + d.p(0)/m.hp(0) )
# alphab = k/h[0]*( -ap(0)*h[1]/(2*hp(0)) - dp(0,kappa)/hp(0) )
# alphac = k/h[0]*( -ap(0)*h[0]/(2*hp(0)) + dp(0,kappa)/hp(0) )


# Right hand side Robin boundary coefficients for matrix equation
robin_aJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-1)/(2*m.hm(m.J-1)) + d.m(m.J-1)/m.hm(m.J-1) )
robin_bJ = lambda a, d, m, k: k/m.h(m.J-1)*( a.m(m.J-1)*m.h(m.J-2)/(2*m.hm(m.J-1)) - d.m(m.J-1)/m.hm(m.J-1) )
# betaa = k/h[J-1]*( am(J-1)*h[J-1]/(2*hm(J-1)) + dm(J-1,kappa)/hm(J-1) )
# betab = k/h[J-1]*( am(J-1)*h[J-2]/(2*hm(J-1)) - dm(J-1,kappa)/hm(J-1) )


# Vector elements for Robin boundary condition
gR_L, gR_R = (0, 0) # left and right Robin boundary values
robin_bcv1 = lambda k, m, gR1: k*gR1/m.h(0)
robin_bcvJ = lambda k, m, gRJ: k*gRJ/m.h(m.J-1)


# # Vector elements for Dirichlet boundary condition
# gD_L, gD_R = (1, 0)     # left and right Dirichlet boundary values
# bL_D = gD_L             # @ left boundary
# bL1_D = ra(1)*gD_L      # @ left boundary + 1
# b1R_D = rc(J-2)*gD_R    # @ right boundary - 1
# bR_D = gD_R             # @ right boundary




# # Delete the work directory (for output files)
# import os, shutil
# working_directory = os.path.join("/tmp", "output")
# if os.path.exists(working_directory):
#     # This is for safey, we will only ever delete folders in the tmp directory
#     dirs = os.path.split(working_directory)
#     if dirs[0] == "/tmp" and len(dirs)>1:
#         shutil.rmtree(working_directory)
#     else:
#         print "Maybe the working directory %s is not safe to delete. Aborting."

# Re-create the working directory
# if not os.path.exists(working_directory):
#     os.makedirs(working_directory)

# # Left are right corner values for Dirichlet BCs
# left_corner = (1,0,0); right_corner = (0,0,1)
# A = advection_diffusion_A_matrix(ra, rb, rc, left_corner, right_corner, J)
# left_corner = (0,0,0); right_corner = (0,0,0)
# M = advection_diffusion_M_matrix(ra, rb, rc, left_corner, right_corner, J)
# b = np.zeros(len(cells))
# b[0] = bL_D
# b[1] = bL1_D
# b[-2] = b1R_D
# b[-1] = bR_D

# # Analytical solution for Dirichlet boundary conditions
# analytical_x = np.concatenate([np.array([0]), cells, np.array([1])])
# analytical_solution = np.concatenate([np.array([1]), (np.exp(a/d) - np.exp(cells*a/d))/(np.exp(a/d)-1), np.array([0]) ])


# Left are right corner values for Robin BCs
left_corner = (1-theta*robin_b1(a,d,mesh,k), -theta*robin_c1(a,d,mesh,k), -theta*ra(1,a,d,mesh,k))
right_corner = (-theta*rc(J-2,a,d,mesh,k), -theta*robin_aJ(a,d,mesh,k),1 - theta*robin_bJ(a,d,mesh,k))
A = advection_diffusion_A_matrix(a, d, mesh, k, ra, rb, rc, left_corner, right_corner, J)

left_corner = (1+(1-theta)*robin_b1(a,d,mesh,k),(1-theta)*robin_c1(a,d,mesh,k),(1-theta)*ra(1,a,d,mesh,k)); 
right_corner = ((1-theta)*rc(J-2,a,d,mesh,k), (1-theta)*robin_aJ(a,d,mesh,k), 1+(1-theta)*robin_bJ(a,d,mesh,k))
M = advection_diffusion_M_matrix(a, d, mesh, k, ra, rb, rc, left_corner, right_corner, J)
b = np.zeros(mesh.J)
b[0] = robin_bcv1(k,mesh,gR_L)
b[-1] = robin_bcvJ(k,mesh,gR_R)

# Source term
b[int(np.median(range(mesh.J)))] = 0.01

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
plt.ylim(0,1.1)
l1.set_data(mesh.cells,w_init)

w = w_init
with writer.saving(fig, "writer_test.mp4", 300):
    
    for i in range(3001):
        b[-1] = w[-1] * -0.1
        b[int(np.median(range(mesh.J)))] = np.random.uniform(0.1,0)
        w = linalg.spsolve(A.tocsc(), M * w + b)
        
        if  i == 0:
            l1.set_data(mesh.cells,w_init)
            writer.grab_frame()
            
        if i %  20 == 0 or i == 0:
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



    