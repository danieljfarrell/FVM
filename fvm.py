from __future__ import division
import numpy as np
import pylab
from scipy.sparse import linalg

faces = np.linspace(-0.5,1,300)   # Position of cell faces
#faces = np.array([0,1,2,3])
#faces = np.array([-2,1,4,8])
cells = 0.5 * (faces[0:-1] + faces[1:]) # Position of cell centroids
J = len(cells)                  # Total number of cells
L = np.max(faces) - np.min(faces)
a = 0.1*np.ones(len(cells))     # Advection velocity
d = 1e-3*np.ones(len(cells))     # Diffusion coefficient
k = 0.001                      # Time step 
theta = 0.5                     # Implicit/Explicit


# Define spacing between cell faces and centroids
h = faces[1:] - faces[0:-1]
hm = lambda i: cells[i] - cells[i-1] # We can use numpy array arguments
hp = lambda i: cells[i+1] - cells[i] # here without issue (vectorisable!)
#print "cell widths", h
#print "backward centroid separations", hm(np.array(range(1,J)))
#print "forward centroid separations", hp(np.array(range(0,J-1)))


# Velocity at cell's left and right hand faces, linear interpolation
am = lambda i: h[i]/(2*hm(i))*a[i-1] + h[i-1]/(2*hm(i))*a[i]
ap = lambda i: h[i+1]/(2*hp(i))*a[i] + h[i]/(2*hp(i))*a[i+1]
#print "velocity at faces using backwards linear interpolation", am(np.array(range(1,J)))
#print "velocity at faces using forwards linear interpolation", ap(np.array(range(0,J-1)))


# Diffusion coefficient at cell's left and right hand faces, cell averages
dm = lambda i: h[i]/(2*hm(i))*d[i-1] + h[i-1]/(2*hm(i))*d[i]
dp = lambda i: h[i+1]/(2*hp(i))*d[i] + h[i]/(2*hp(i))*d[i+1]
#print "diffusion coefficient at faces using backwards linear interpolation", dm(np.array(range(1,J)))
#print "diffusion coefficient at faces using forwards linear interpolation", dp(np.array(range(0,J-1)))


# Interior coefficients for matrix equation
ra = lambda i: k/h[i]*(am(i)*h[i]/(2*hm(i)) + dm(i)/hm(i))
rb = lambda i: k/h[i]*(am(i)*h[i-1]/(2*hm(i)) - ap(i)*h[i+1]/(2*hp(i)) - dm(i)/hm(i) - dp(i)/hp(i))
rc = lambda i: k/h[i]*(-ap(i)*h[i]/(2*hp(i)) + dp(i)/hp(i))
#print "'r' coefficients evaluates at j=2:"
#print "\tra", ra(1) 
#print "\trb", rb(1) 
#print "\trc", rc(1) 

# Left hand side Robin boundary coefficients for matrix equation
alphab = k/h[0]*( -ap(0)*h[1]/(2*hp(0)) - dp(0)/hp(0) )
alphac = k/h[0]*( -ap(0)*h[0]/(2*hp(0)) + dp(0)/hp(0) )
#print "'alpha' boundary terms at x=xL"
#print "\talphab", alphab
#print "\talphac", alphac

# Right hand side Robin boundary coefficients for matrix equation
betaa = k/h[-1]*( am(-1)*h[-1]/(2*hm(-1)) + dm(-1)/hm(-1) )
betab = k/h[-1]*( am(-1)*h[-2]/(2*hm(-1)) - dm(-1)/hm(-1) )
#print "'alpha' boundary terms at x=xL"
#print "\tbetaa", betaa
#print "\tbetab", betab

# Now write in matrix from, theta-method has been applied

# Left hand side matrix
from scipy import sparse
padding = np.array([0])
inx = np.array(range(1,J-1))
central = np.concatenate([np.array([1-theta*alphab]), 1-theta*rb(inx), np.array([1-theta*betab])])
inx = np.array(range(1,J-1))
lower = np.concatenate([-theta*ra(inx), np.array([-theta*betaa]), padding])
inx = np.array(range(1,J-1))
upper = np.concatenate([padding, np.array([-theta*alphac]), -theta*rc(inx)])
diags = [-1,0,1]
A = sparse.spdiags([lower, central, upper], diags, J, J)
#print "central", central
#print "lower", lower
#print "upper", upper
#print A.todense()


# Right hand side matrix 
padding = np.array([0])
inx = np.array(range(1,J-1))
central = np.concatenate([np.array([1+(1-theta)*alphab]), 1+(1-theta)*rb(inx), np.array([1+(1-theta)*betab])])
inx = np.array(range(1,J-1))
lower = np.concatenate([(1-theta)*ra(inx), np.array([(1-theta)*betaa]), padding])
inx = np.array(range(1,J-1))
upper = np.concatenate([padding, np.array([(1-theta)*alphac]), (1-theta)*rc(inx)])
diags = [-1,0,1]
M = sparse.spdiags([lower, central, upper], diags, J, J)
#print "central", central
#print "lower", lower
#print "upper", upper
#print M.todense()


# Check Peclet number and CFL condition
mu = a * h / d
CFL = a * k / h
print "Peclet number", np.average(mu)
print "CFL condition", np.average(CFL)
#pylab.plot(cells, mu, label="Peclet number")
#pylab.plot(cells, a*k/h, label="CFL condition")
#pylab.ylim((-1,None))
#pylab.legend()
#pylab.show()


# Initial conditions
gaussian = lambda z, height, position, hwhm: height * np.exp(-np.log(2) * ((z - position)/hwhm)**2)
w_init = gaussian(cells, 1, 0, 0.1)
#w_init = np.where(w_init < 1e-5, 0, w_init)
#pylab.plot(cells, w_init)
#pylab.show()

area_0 = np.trapz(w_init, x=cells)
sum_0 = np.sum(w_init) * L/J
print area_0, sum_0

import os
directory = os.path.join(os.path.dirname(os.path.realpath(__file__)), "output")

# Time sweep
w = w_init
for i in range(30000):
    w = linalg.spsolve(A.tocsc(), M * w)
    #w = linalg.gmres(A.tocsc(), M*w)[0]
    #w = linalg.lgmres(A.tocsc(), M*w)[0]
    #w = linalg.minres(A.tocsc(), M*w)[0]
    #print len(w), w
    if i %  1000 == 0:
        pylab.ylim((0,10))
        pylab.fill_between(cells, 0, w)
        pylab.plot(cells, w)
        area = np.trapz(w, x=cells)
        summation = np.sum(w) * L/J
        print area, summation
        pylab.title("Area: (numpy.trapz) %g; (simple summation) %g" % (area, summation))
        #pylab.plot(cells, a*h/d)
        #pylab.draw()
        pylab.savefig(os.path.join(directory, "%d_solution.png" % i), dpi=72)
        pylab.cla()

def covert_images_to_animated_gif(target_dir, filename="animation2.gif"):
    from images2gif import writeGif
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
    writeGif(filename, images, duration=0.1)
    
covert_images_to_animated_gif(directory, filename="solution.gif")


