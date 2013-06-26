from __future__ import division
import collections
import numpy as np
import pylab
from scipy import sparse
from scipy.sparse import linalg
np.random.seed(seed=1)

# Supporting functions
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
        writeGif(filename, images, duration=0.2)
    




left_face = lambda i, f, h, hm: h[i]/(2*hm(i))*f[i-1] + h[i-1]/(2*hm(i))*f[i] # Interpolation from cell to face values
right_face = lambda i, f, h, hp: h[i+1]/(2*hp(i))*f[i] + h[i]/(2*hp(i))*f[i+1]
    
faces = np.linspace(-0.5,1,50)   # Position of cell faces
#faces = np.concatenate((np.array([-0.5]), np.sort(np.random.uniform(-0.5, 1, 50)), np.array([1])))
cells = 0.5 * (faces[0:-1] + faces[1:]) # Position of cell centroids
J = len(cells)                  # Total number of cells
L = np.max(faces) - np.min(faces)

a = 0.1*np.ones(len(cells))     # Advection velocity
#a = np.where( (cells > 0)*(cells<0.5), a, 0)
d = 0.04*np.ones(len(cells))    # Diffusion coefficient
k = 0.01                       # Time step 
theta = 0.5                     # Implicit/Explicit


# Define mesh face and centroid spacings
h = faces[1:] - faces[0:-1]
def hm(i):
    if not check_index_within_bounds(i,1,J-1):
        raise ValueError("hm index runs out of bounds")
    return cells[i] - cells[i-1]

def hp(i):
    if not check_index_within_bounds(i,0,J-2):
        raise ValueError("hp index runs out of bounds")
    return cells[i+1] - cells[i]





# Velocity on the cell faces
am = lambda i: left_face(i, a, h, hm)
ap = lambda i: right_face(i, a, h, hp)


# Diffusion coefficient on the cell faces
dm = lambda i, kappa: h[i]/(2*hm(i))*d[i-1] + h[i-1]/(2*hm(i))*d[i] + 0.5 * kappa[i] * hm(i) * am(i)
dp = lambda i, kappa: h[i+1]/(2*hp(i))*d[i] + h[i]/(2*hp(i))*d[i+1] + 0.5 * kappa[i] * hp(i) * ap(i)

# Check Peclet number and CFL condition
mu = a * h / d
CFL = a * k / h
print "Peclet number", np.average(mu), np.max(mu)
print "CFL condition", np.average(CFL), np.max(CFL)
kappa = np.sign(a) * 0

#kappa = (np.exp(mu) + 1)/(np.exp(mu) - 1) - 2/mu
kappa = np.zeros(len(mu)) # Force central difference


# Interior coefficients for matrix equation
ra = lambda i: k/h[i]*(am(i)*h[i]/(2*hm(i)) + dm(i,kappa)/hm(i))
rb = lambda i: k/h[i]*(am(i)*h[i-1]/(2*hm(i)) - ap(i)*h[i+1]/(2*hp(i)) - dm(i,kappa)/hm(i) - dp(i,kappa)/hp(i))
rc = lambda i: k/h[i]*(-ap(i)*h[i]/(2*hp(i)) + dp(i,kappa)/hp(i))


# Left hand side Robin boundary coefficients for matrix equation
alphab = k/h[0]*( -ap(0)*h[1]/(2*hp(0)) - dp(0,kappa)/hp(0) )
alphac = k/h[0]*( -ap(0)*h[0]/(2*hp(0)) + dp(0,kappa)/hp(0) )


# Right hand side Robin boundary coefficients for matrix equation
betaa = k/h[J-1]*( am(J-1)*h[J-1]/(2*hm(J-1)) + dm(J-1,kappa)/hm(J-1) )
betab = k/h[J-1]*( am(J-1)*h[J-2]/(2*hm(J-1)) - dm(J-1,kappa)/hm(J-1) )


# Vector elements for Robin boundary condition
gR_L, gR_R = (0, 0) # left and right Robin boundary values
bL_R = k*gR_L/h[0]
bR_R = k*gR_R/h[J-1]


# Vector elements for Dirichlet boundary condition
gD_L, gD_R = (1, 0)     # left and right Dirichlet boundary values
bL_D = gD_L             # @ left boundary
bL1_D = ra(1)*gD_L      # @ left boundary + 1
b1R_D = rc(J-2)*gD_R    # @ right boundary - 1
bR_D = gD_R             # @ right boundary


# # Dirichlet boundary conditions, left hand side
# gD = 0
# alphab_D = k/h[0]*( -ap(0)*h[1]/(2*hp(0)) - 2*dm(1,kappa)/h[0] - dp(0,kappa)/hp(0) ) # NB dm(0) will not evaluate we need to interpolate to this value
# alphac_D = k/h[0]*( -ap(0)*h[0]/(2*hp(0)) + dp(0,kappa)/hp(0) )
# bL_D = k*gD/h[0] * (am(1) + 2*dm(1,kappa)/h[0] ) # NB am(0) and dm(0) will not evaluate we need to interpolate to this value


# # Dirichlet boundary conditions, right hand side
# gD = 1
# j = J-1
# betab_D = k/h[j]*( -am(j)*h[j-1]/(2*hm(j)) + 2*dp(j-1,kappa)/h[j] + dm(j,kappa)/hm(j) ) # NB dp(j) (dp @ right boundary), will not evaluate we need to interpolate to this value
# betaa_D = k/h[j]*( -am(j)*h[j]/(2*hm(j)) - dm(j,kappa)/hm(j) )
# bR_D = k*gD/h[j] * (ap(j-1) - 2*dp(j-1,kappa)/h[j] ) # NB ap(j) and dp(j) (ap and dp @ right boundary) will not evaluate we need to interpolate to this value


# Pick which boundary conditions term you want and choose terms for the boundary condition vector

# # Robin boundary conditions
# alphab, alphac = alphab_R, alphac_R
# betaa, betab = betaa_R, betab_R
# bL, bR = bL_R, bR_R

# Left and right hand side matrices using theta-method
# def A_matrix_Robin_BCs(ra, rb, rc, alphab, alphac, betaa, betab, J):
#     padding = np.array([0])
#     inx = np.array(range(1,J-1))
#     upper = np.concatenate([padding, np.array([-theta*alphac]), -theta*rc(inx)])
#     central = np.concatenate([np.array([1-theta*alphab]), 1-theta*rb(inx), np.array([1-theta*betab])])
#     lower = np.concatenate([-theta*ra(inx), np.array([-theta*betaa]), padding])
#     return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

# def M_matrix_Robin_BCs(ra, rb, rc, alphab, alphac, betaa, betab, J):
#     padding = np.array([0])
#     inx = np.array(range(1,J-1))
#     upper = np.concatenate([padding, np.array([(1-theta)*alphac]), (1-theta)*rc(inx)])
#     central = np.concatenate([np.array([1+(1-theta)*alphab]), 1+(1-theta)*rb(inx), np.array([1+(1-theta)*betab])])
#     lower = np.concatenate([(1-theta)*ra(inx), np.array([(1-theta)*betaa]), padding])
#     return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

# def A_matrix_Dirichlet_BCs(ra, rb, rc, J):
#     padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
#     zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
#     one  = np.array([1])    #
#     inx = np.array(range(1,J-2))
#     upper = np.concatenate([padding, zero, -theta*rc(inx), zero])
#     inx = np.array(range(1,J-1))
#     central = np.concatenate([one, 1-theta*rb(inx), one])
#     inx = np.array(range(2,J-1))
#     lower = np.concatenate([zero, -theta*ra(inx), zero, padding])
#     return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

# def M_matrix_Dirichlet_BCs(ra, rb, rc, J):
#     padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
#     zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
#     one  = np.array([1])    #
#     inx = np.array(range(1,J-2))
#     upper = np.concatenate([padding, zero, (1-theta)*rc(inx), zero])
#     inx = np.array(range(1,J-1))
#     central = np.concatenate([zero, 1+(1-theta)*rb(inx), zero])
#     inx = np.array(range(2,J-1))
#     lower = np.concatenate([zero, (1-theta)*ra(inx), zero, padding])
#     return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

def advection_diffusion_A_matrix(ra, rb, rc, left_corner, right_corner, J):
    padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
    zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
    one  = np.array([1])    #
    inx = np.array(range(1,J-2))
    b1, c1, a2 = [np.array([item]) for item in left_corner]
    c1J, aJ, bJ = [np.array([item]) for item in right_corner]
    upper = np.concatenate([padding, c1, -theta*rc(inx), c1J])
    inx = np.array(range(1,J-1))
    central = np.concatenate([b1, 1-theta*rb(inx), bJ])
    inx = np.array(range(2,J-1))
    lower = np.concatenate([a2, -theta*ra(inx), aJ, padding])
    return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)

def advection_diffusion_M_matrix(ra, rb, rc, left_corner, right_corner, J):
    padding = np.array([0]) # A element which is pushed off the edge of the matrix by the spdiags function
    zero = padding          # Yes, its the same. But this element is included in the matrix (semantic difference).
    one  = np.array([1])    #
    inx = np.array(range(1,J-2))
    b1, c1, a2 = [np.array([item]) for item in left_corner]
    c1J, aJ, bJ = [np.array([item]) for item in right_corner]
    upper = np.concatenate([padding, c1, (1-theta)*rc(inx), c1J])
    inx = np.array(range(1,J-1))
    central = np.concatenate([b1, 1+(1-theta)*rb(inx),bJ])
    inx = np.array(range(2,J-1))
    lower = np.concatenate([a2, (1-theta)*ra(inx), aJ, padding])
    return sparse.spdiags([lower, central, upper], [-1,0,1], J, J)




# Initial conditions
gaussian = lambda z, height, position, hwhm: height * np.exp(-np.log(2) * ((z - position)/hwhm)**2)
H = lambda z: 0.5 * (1 - np.sign(z))
#w_init = gaussian(cells, 1, 0, 0.1)
#w_init[0] = 0; w_init[-1] = 0
w_init = H(cells)
#w_init = np.where(w_init < 1e-5, 0, w_init)
#pylab.plot(cells, w_init)
#pylab.show()


# Delete the work directory (for output files)
import os, shutil
working_directory = os.path.join("/tmp", "output")
if os.path.exists(working_directory):
    # This is for safey, we will only ever delete folders in the tmp directory
    dirs = os.path.split(working_directory)
    if dirs[0] == "/tmp" and len(dirs)>1:
        shutil.rmtree(working_directory)
    else:
        print "Maybe the working directory %s is not safe to delete. Aborting."

# Remake the working directory
if not os.path.exists(working_directory):
    os.makedirs(working_directory)


# Time sweep
w = w_init

#b[0] = bL; b[-1] = bR
# A = A_matrix_Robin_BCs(ra,rb,rc,alphab,alphac,betaa,betab,J)
# M = M_matrix_Robin_BCs(ra,rb,rc,alphab,alphac,betaa,betab,J)
#A = A_matrix_Dirichlet_BCs(ra,rb,rc,J)

# Left are right corner values for Dirichlet BCs
left_corner = (1,0,0); right_corner = (0,0,1)
A = advection_diffusion_A_matrix(ra, rb, rc, left_corner, right_corner, J)
left_corner = (0,0,0); right_corner = (0,0,0)
M = advection_diffusion_M_matrix(ra, rb, rc, left_corner, right_corner, J)
b = np.zeros(len(w))
b[0] = bL_D
b[1] = bL1_D
b[-2] = b1R_D
b[-1] = bR_D

# Left are right corner values for Robin BCs
left_corner = (1-theta*alphab, -theta*alphac, -theta*ra(1)); 
right_corner = (-theta*rc(J-2), -theta*betaa,1 - theta*betab)
A = advection_diffusion_A_matrix(ra, rb, rc, left_corner, right_corner, J)
left_corner = (1+(1-theta)*alphab,(1-theta)*alphac,(1-theta)*ra(1)); 
right_corner = ((1-theta)*rc(J-2), (1-theta)*betaa, 1+(1-theta)*betab)
M = advection_diffusion_M_matrix(ra, rb, rc, left_corner, right_corner, J)
b = np.zeros(len(w))
b[0] = bL_R
b[-1] = bR_R

#analytical_x = np.concatenate([np.array([0]), cells, np.array([1])])
#analytical_solution = np.concatenate([np.array([1]), (np.exp(a/d) - np.exp(cells*a/d))/(np.exp(a/d)-1), np.array([0]) ])


# import pylab
# cmap = pylab.cm.Blues
# pylab.matshow(A.todense(),cmap=cmap)
# pylab.show()
# pylab.cla()
# pylab.matshow(M.todense(),cmap=cmap)
# pylab.show()
# pylab.cla()
# pylab.matshow(np.matrix(b), cmap=cmap)
# pylab.show()
# pylab.cla()
#exit()

for i in range(4000):
    w = linalg.spsolve(A.tocsc(), M * w + b)
    if i %  10 == 0 or i == 0:
        pylab.ylim((0,2))
        pylab.bar(faces[:-1], w, width=(faces[1:]-faces[:-1]), edgecolor='none')
        pylab.plot(cells, w, "k", lw=3)
        #pylab.plot(analytical_x, analytical_solution, "r-", lw=2)
        area = np.sum(w * h)
        print "#%d; area:" % (i,), area
        pylab.title("Area: %g" % (area))
        pylab.savefig(os.path.join(working_directory, "%d_solution.png" % i), dpi=72)
        pylab.cla()

# Write the animated gif to the script directory
pwd = os.path.join(os.path.dirname(os.path.realpath(__file__)))
covert_images_to_animated_gif(working_directory, filename=os.path.join(pwd, "solution.gif") )