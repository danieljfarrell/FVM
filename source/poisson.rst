The finite volume method applied to the Poisson equation
--------------------------------------------------------

The advection-diffusion equation used in electrostatics is of the form,

.. math::
	0 = (\epsilon(x) \phi_x)_x + \rho(x,t)

where :math:`\phi` is the electric potential and the :math:`x` subscripts mean a spatial partial derivative. Note that this equation is reaction-diffusion equation with the diffusion flux :math:`\mathcal{F} = -\epsilon(x) \phi_x`. As with the advection-diffusion equation, this expression holds over the whole domain, therefore it must also hold over a subdomain :math:`\Omega_j=(x_{j-1/2}, x_{j+1/2})`. Integrating over the subdomain,
	
.. math::
	0 = \int_{x_{j-1/2}}^{x_{j+1/2}} (\mathcal{-F})_x~dx + \int_{x_{j-1/2}}^{x_{j+1/2}} \rho(x,t)~dx

The cell averaged values are therefore,

.. math::
	0 =  -\frac{\mathcal{F}_{j+1/2}}{h_j} + \frac{\mathcal{F}_{j-1/2}}{h_{j}} + \bar{\rho}_j

To determine the flux at the cell faces we will take the mean of the fluxes in cells adjacent to the interfaces,

.. math::
	\mathcal{F}_{j+\frac{1}{2}} = - d_{j+\frac{1}{2}} \frac{\phi_{j+1}-\phi_j}{h_{+}}

.. math::
	\mathcal{F}_{j-\frac{1}{2}} = - d_{j-\frac{1}{2}} \frac{\phi_{j}-\phi_{j-1}}{h_{-}}
  
Substituting the definition of the fluxes into the integrated equation and factoring for the electrical potential terms yields,

.. math::
	0 = \frac{\phi_{j-1}}{h_j} \left( \frac{\epsilon(x_{j-\frac{1}{2}})}{h_{-}}  \right) - \frac{\phi_{j}}{h_j} \left( \frac{\epsilon(x_{j-\frac{1}{2}})}{h_{-}} + \frac{\epsilon(x_{j+\frac{1}{2}})}{h_{+}}  \right) + \frac{\phi_{j+1}}{h_j} \left( \frac{\epsilon(x_{j+\frac{1}{2}})}{h_{+}} \right) + \bar{\rho}_j

Naturally resolved Neumann boundary conditions
----------------------------------------------

Neumann boundary conditions applied to the finite volume form of the Poisson equation are naturally resolved; that is to say that ghost cell or extrapolation approaches are not needed. Consider the integral form of the Poisson equation at the last cell of the domain on the right hand side,

.. math::
	0 = \frac{1}{h_J}\left( \epsilon(x_{R}) \phi_x\bigg|_{x=x_R} - \epsilon(x_{J-1/2}) \phi_x\bigg|_{J-1/2} \right) + \bar{\rho}_J

The Neumann boundary condition defines the derivative of :math:`\phi` at the the right hand boundary,

.. math::
	\phi_x\bigg|_{x=x_R} = \sigma_R

This boundary condition appears naturally as a term in the integral equation. Therefore the boundary condition can be satisfied by substitution into the integral form of the equation at the right hand boundary,

.. math::
	\frac{1}{h_J}\left( \epsilon(x_{R}) \sigma_R - \epsilon(x_{J-1/2}) \phi_x\bigg|_{J-1/2} \right) + \bar{\rho}_J

Dirichlet boundary conditions
-----------------------------

Fixed value, or Dirichlet, boundary conditions are not naturally resolved and can only be applied to the discretised form of the equation. We need to apply a Dirichlet condition to the left hand side, such that full equation for the cell :math:`j=1` reads,

.. math::
	\phi(x_L) = \omega_L

Clearly, this is trivial as all that is required is fix the first row of the matrix equation to be the boundary condition.

Poisson equation in matrix form
-------------------------------

Defining the coefficients,

.. math::
	\ell_a & = \frac{1}{h_j}\frac{\epsilon(x_{j-1/2})}{h_{-}} \\
	\ell_b & = -\frac{1}{h_j}\left( \frac{\epsilon_{j-1/2}}{h_{-}} + \frac{\epsilon_{j+1/2}}{h_{+}} \right) \\
	\ell_c & = \frac{1}{h_j}\frac{\epsilon(x_{j+1/2})}{h+} \\

With Dirichlet boundary condition on the left hand side and a Neumann boundary condition on the right hand side the Poisson equation can be written in matrix form,

.. math::
	\begin{align} 
	\begin{pmatrix}
	1			& 0    		&        		&       		& 0   			\\
	\ell_a(2)	& \ell_b(2)	& \ell_c(2)    	&       		&     			\\
	    		& \ddots 	& \ddots 		& \ddots		&     			\\
	    		&        	&  \ell_a(J-1) 	& \ell_b(J-1)   & \ell_c(J-1)	\\
	 0  		&        	&        		& \ell_a(J)   	& -\ell_a(J)	
	\end{pmatrix}
	\begin{pmatrix}
	    \phi_1 \\
	    \phi_2 \\
	    \vdots \\
	    \phi_{J-1} \\
	    \phi_{J} \\
	\end{pmatrix}
	+
	\begin{pmatrix}
	    -\omega_L \\
	    \bar{\rho}_2 \\
	    \vdots \\
	    \bar{\rho}_{J-1} \\
	    1/h_j\epsilon(x_J)\sigma_R + \bar{\rho}_{J} \\
	\end{pmatrix}
	= 0
	\end{align}
	



