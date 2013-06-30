Solver implementation in Python
-------------------------------

Here we discuss how to implement a solver for the advection-diffusion equation in Python. The notes will consider how to design a solver which minimises code complexity and maximise readability.

Overview
********

At a high-level usage of the code looks like the following,

.. code-block:: python
	:emphasize-lines: 3,6,7,11
	
	# Define a mesh
	faces = np.linspace(-0.5, 1, 100)
	mesh =  Mesh(faces)
	
	# Define coefficients
	a = CellVariable(0.01, mesh=mesh) # Advection velocity
	d = CellVariable(1e-3, mesh=mesh) # Diffusion coefficient
	
	# Make a 'model' and apply boundary conditions
	k = 1 # Time step
	model = Model(faces, a, d, k)
	model.set_boundary_conditions(left_value=1., right_value=0.)
	
	# Ask for the linear system... 
	A = model.A_matrix()
	M = model.M_matrix()
	b = model.b_vector()
	
	# Solve...

There are a number of classes (highlighted above) which abstract away details of dealing with finite-volume equations:

 1. :code:`Mesh`
 2. :code:`CellVariable`
 3. :code:`Model`

The the :code:`Mesh` and :code:`CellVariable` classes have been inspired by the approach of Fipy_.

.. _Fipy: http://www.ctcms.nist.gov/fipy/

In the following we will highlight the main features of each class and point out how they are useful. The classes do not attempt at doing "*too much*", they simply aid in the creation of the linear system.

The :code:`Mesh` class
**********************

:code:`Mesh` objects are initialised with a list of faces locations, which can be non-uniformly distributed if desired. A mesh is completely defined by locations of cell **faces**. Some useful methods,

.. code:: python
	
	def h(self, i):
		...

:code:`h()` returns the cell width for cell with index :code:`i`. This is function is vectorisable by passing an array of the required indices but it *does not* accept "fancy indexing". The reason being, that it is very easy to make mistakes with subscript indexing, the goal here is to make the user be explicit when requesting elements. Note the :code:`self.cell_widths` instance variable returns the numpy array of cell widths is a second way of accessing this data.

.. code:: python

	def hm(self, i):
	...

:code:`hm()` returns the distance between cell centroids for the cell at index :code:`i` and :code:`i-1`, that is in the *backwards* or **minus** direction.

.. code:: python

	def hp(self, i):
	...

Similarly :code:`hp()` returns the distance between cell centroids for the cell at index :code:`i` and :code:`i+1`, that is in the *forwards* or **plus** direction.

In addition to the above method, the class contains the instance variables, :code:`self.cells` which returns an array of cell centroid locations, :code:`self.J` which contains the number of cells, and also a copy of the face locations (an array) via :code:`self.faces`.

The :code:`CellVariable` class
******************************

The goal of the :code:`CellVariable` class is to provide a elegant way of automatically interpolating between the cell value and the face value. The class holes values which correspond to the **cell average**. Internally, this class is a subclass of :code:`numpy.ndarray` so it is a fully functioning numpy array. It has a new constructor and additional method which return interpolated values at the cell surfaces.

A :code:`CellVariable` is initialised with a value for the cell average (this can be a constant or an array-like quantity) and the :code:`Mesh` on which the cell variable is defined. My coupling the cell variable with the mesh the class can perform interpolation between the cell and face values using the methods,

.. code:: python

	def p(self, i):
	...
	
	def m(self, i)
	...

Again :code:`self.p(i)` stands for the *plus* direction and :code:`self.m(i)` stands for the *minus* direction, as such they return values at the right and left face of the cell. The mesh can be returned via the instance variable :code:`cell_variable.mesh`.


The :code:`Model` class
***********************

The model class is where the creating of the matrices occurs and where boundary conditions can be applied to the problem. For these reasons the class is fairly complicated.

There are method which return different element of the final matrix. The interior elements are fairly homogenous, the only real difference is where there are spatially varying coefficient of cell widths. For this reason the the method :code:`_interior_elements()` returns elements which correspond to the lower, central and upper diagonals for a specific index. For example, to calculate the interior matrix elements for mesh point :code:`i=4` one would do the following,

.. code:: python

	model = Model(...)
	ra, rb, rc = model._interior_functions()
	"index is i=4"
	ra(4, model.a, model.d, model.mesh, model.k) # lower diagonal function
	ra(4, model.a, model.d, model.mesh, model.k) # central diagonal function
	ra(4, model.a, model.d, model.mesh, model.k) # upper diagonal function

The function names here correspond to the matrix element in the previous section.

Note that the function is prefixed with an underscore this is because a 'users' should have no need to call this method. It is called internally when constructing the finite-volume matrices. However, an 'author' does need to provide the correct matrix element with this function.

The methods,

.. code:: python

	def _robin_boundary_condition_elements_left(self):
		...
		
 	def _robin_boundary_condition_elements_right(self):
		...
		
 	def _dirichlet_boundary_condition_elements_left(self):
		...
		
	def _dirichlet_boundary_condition_elements_right(self):
		...

play a similar role. However the return a list of index-value pairs :code:`([(1,1), a11], [(i,2), b12] ...)` rather than returning functions. The functions return the value of element which need to change (with respect to the interior values) in order include boundary conditions. The index-value pair facilitates automatic insertion of the values into the correct matrix element. As we will see later, rather than hard coding the position of the various element if the index and value are specified it makes the destination of the element unambiguous. It also allows the value of the matrix element to be defined at the same point in the code as the location. This is beneficial for providing context and should reduce bugs and complexity.
 
Boundary conditions modify terms in the :math:`\boldsymbol{A}` and :math:`\boldsymbol{M}` matrices by they also require that a vector be added to the equations. The form of the linear system being solved is,

.. math::
	\boldsymbol{A} \cdot w^{n+1} = \boldsymbol{M} \cdot w^n = b

where :math:`b` is a vector contains the boundary conditions values (and also values of the source term should it exist). The elements of :math:`b` are returned from the following methods, 

.. code:: python

	def _robin_boundary_condition_vector_elements_left(self):
		...
		
 	def _robin_boundary_condition_vector_elements_right(self):
		...
		
 	def _dirichlet_boundary_condition_vector_elements_left(self):
		...
		
	def _dirichlet_boundary_condition_vector_elements_right(self):
		...
 
Again, these method should return *index-values* pairs, but because the are element of a vector the index is simply a number, not a tuple as with the matrix elements.

The :code:`Model` class also include some convenience function for checking the value of the Peclet number and the CFL conditions which can be called via,

.. code:: python

	def peclet_number(self):
		return self.a * self.mesh.cell_widths / self.d
	
   	def CFL_condition(self):
		return self.a * self.k / self.mesh.cell_widths
		

The method which are intended for the user to actually call when constructing the linear system are,


.. code:: python

   def A_martrix(self):
   		...
    
   def M_martrix(self):
   		...

   def b_vector(self):
   		...

Which simply return the matrices and vector of the linear system.

Finally, when initialising a :code:`Model` object two important keyword arguments can be passed, they are, :code:`theta` and :code:`, discretisation`. The value of :code:`theta` controls the time-integration method (setting :code:`theta=0.5` achieved a Crank-Nicolson trapezoidal integration in time), and the value of :code:`discretisation` can be one of the following: :code:`'upwind', 'central', 'exponential'`. The :code:`upwind` option uses the classic *first order upwind* discretisation, :code:`central` uses *second-order central* and setting to :code:`exponential` uses an adaptive scheme which will use weight between the central and upwind scheme depending on the local value of the Peclet number. This is the classic 'exponential fitting' or 'Scharfetter-Gummel' discretisation. **N.B.** Scharfetter-Gummel also refers to a method of solving the advection-diffusion equation is a non-coupled manner, this is not the case here where it only refers to the the discretisation method.
