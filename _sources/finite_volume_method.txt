The finite volume method
------------------------

The advection-diffusion equation with non-linear source term, in flux form, is given by,

.. math::
	u_t = (\mathcal{-F})_x + s(t,x,u)

where :math:`\mathcal{F} = au - du_x`. This expression holds over the whole domain, therefore it must also hold over a subdomain :math:`\Omega_j=(x_{j-1/2}, x_{j+1/2})`. Integrating over the subdomain, or cell, gives,
	
.. math::
	\int_{x_{j-1/2}}^{x_{j+1/2}} u_t~dx = \int_{x_{j-1/2}}^{x_{j+1/2}} (\mathcal{-F})_x~dx + \int_{x_{j-1/2}}^{x_{j+1/2}} s(x,t,u)~dx

If we now use :math:`w` to represent the cell averages then we don't have to perform the actual integrations,

.. math::
	w_j^{\prime} =  -\frac{\mathcal{F}_{j+1/2}}{h_j} + \frac{\mathcal{F}_{j-1/2}}{h_{j}} + \bar{s}_j

Using linear interpolation (see the next section) to determine the fluxes at the cell faces,

.. math::
	\mathcal{F}_{j+\frac{1}{2}} = a_{j+\frac{1}{2}}\left( \frac{h_{j+1}}{2h_{+}} w_j + \frac{h_j}{2h_{+}} w_{j+1} \right) - d_{j+\frac{1}{2}} \frac{w_j - w_{j+1}}{h_{+}}

.. math::
	\mathcal{F}_{j-\frac{1}{2}} = a_{j-\frac{1}{2}}\left( \frac{h_{j}}{2h_{-}} w_{j-1} + \frac{h_{j-1}}{2h_{-}} w_{j} \right) - d_{j-\frac{1}{2}} \frac{w_{j-1} - w_{j}}{h_{-}}

The meaning of the "h" terms is shown in the figure below,

.. figure:: img/cells_FVM.png
   :scale: 50 %
   :alt: Finite volume cells with important distances and positions labelled.
   :align: center

   Finite volume cells with important distances and positions labelled.
  
Substituting the definition of the fluxes into the above equation yields,

.. math::
	w_j^{\prime} = \frac{w_{j-1}}{h_j} \left( \frac{a(x_{j-\frac{1}{2}}) h_j}{2h_{-}}  - \frac{d(x_{j-\frac{1}{2}})}{h_{-}}\right) + \frac{w_j}{h_j}\left( \frac{a(x_{j-\frac{1}{2}})h_{j-1}}{2h_{-}} - \frac{a(x_{j+\frac{1}{2}})h_{j+1}}{2h_{+}} + 	\frac{d(x_{j-\frac{1}{2}})}{h_{-}} + \frac{d(x_{j+\frac{1}{2}})}{h_{+}}  \right) + \frac{w_{j+1}}{h_j} \left( \frac{-a(x_{j+\frac{1}{2}})h_j}{2h_{+}} - \frac{d(x_{j+\frac{1}{2}})}{h_{+}} \right)

This equation is the *semi-discrete* form.

Adaptive upwinding & exponential fitting
****************************************

Hundsdorfer (on pg. 266) notes that a adaptive upwinding scheme can be introduced by altering the diffusion coefficient,

.. math::
	d(x_{j\pm\frac{1}{2}}) \rightarrow d(x_{j\pm\frac{1}{2}}) + \frac{1}{2}\kappa_{\pm\frac{1}{2}} h_{\pm\frac{1}{2}} a_{j\pm\frac{1}{2}}


.. math::
	\kappa_{\pm}=\text{sgn}(a(x_{j\pm\frac{1}{2}})). 

When :math:`\kappa=0` the discretisation reduces to central different and when :math:`\kappa=\pm1` the discretisation for the advection term becomes an upwind scheme. A generalisation is possible by choosing a function for :math:`\kappa` which automatically adjusts the discretisation method depending on the local value of Peclet's number, :math:`\mu=ah/d`. Exponential fitting has shown to be a robust scheme for steady-state approaches. Exponential fitting can be introduced by,

.. math::
	 \kappa = \frac{e^{\mu}+1}{e^{\mu}-1} - \frac{2}{\mu}

An approximation to the above, which does not require the evaluation of exponentials is commonly used,

.. math::
	\kappa = \begin{cases}
	\text{max}(0, 1-2/\mu) & \text{when}~ \mu>0 \\	
	\text{min}(0, -1-2/\mu) & \text{when}~ \mu<0
	\end{cases}

.. figure:: img/Exponential_fitting_and_approximation.png
   :scale: 10 %
   :alt: Exponential fitting and approximation
   :align: center


This is a very elegant want to introduce adaptive upwinding or exponential fitting because it can be brought into the discretisation in an ad hoc manner. Note that for all of the adaptive expressions :math:`\kappa\rightarrow\pm1` as :math:`\mu\rightarrow\pm\infty`, and :math:`\kappa\rightarrow 0` as :math:`\mu\rightarrow 0`. This means the discretisation will automatically weight in the favour of the upwind scheme in locations where advection dominates. Conversely, where diffusion dominates the weighting factor shifts in favour of a central difference scheme.

Naturally resolved Robin boundary conditions
********************************************

Within the finite volume method Robin boundary conditions are naturally resolved. This means that there is no need for interpolation or ghost point substitution (although these approaches remain possible) to include the boundary conditions because the flux at the boundary appears naturally in the semi-discretised equation. For example consider the semi-discretised equation evaluated at the boundary cell :math:`\Omega_1`,

.. figure:: img/boundary_cell_FVM.png
   :scale: 50 %
   :alt: Finite volume boundary cell at the left hand side.
   :align: center

   Finite volume boundary cell at the left hand side.

.. math::
	w_1^{\prime} =  -\frac{\mathcal{F}_{3/2}}{h_1} + \frac{\mathcal{F}_{1/2}}{h_{1}} + \bar{s}_1

The Robin boundary condition specifies the flux at :math:`x_{1/2}`,

.. math::
	\mathcal{F}_{1/2} = g_{R}(x_{1/2})

Therefore the boundary condition can be incorporated without invoking any information regarding the ghost cell,

.. math::
	w_1^{\prime} = \frac{w_1}{h_1}\left( \frac{-a(x_{3/2})h_{2}}{2h_{+}} + \frac{d(x_{3/2})}{h_{+}} \right) + \frac{w_{2}}{h_1} \left( \frac{-a(x_{3/2}) h_1}{2h_{+}} - \frac{d(x_{3/2})}{h_{+}} \right) + \frac{g_{R}(x_L)}{h_1}

Similarly applying the same procedure to the :math:`\Omega_J` cell at the right hand side boundary,

.. figure:: img/boundary_cell_FVM_rhs.png
   :scale: 50 %
   :alt: Finite volume boundary cell at the right hand side.
   :align: center

   Finite volume boundary cell at the right hand side.

.. math::
	w_J^{\prime} =  -\frac{\mathcal{F}_{J+1/2}}{h_J} + \frac{\mathcal{F}_{J-1/2}}{h_J} + \bar{s}_J

The Robin boundary condition at the right hand side is,

.. math::
	g_{R}(x_R) = \mathcal{F}_{J+1/2}

Therefore the naturally resolved boundary condition on the right hand side becomes,

.. math::
	w_J^{\prime} = \frac{w_{J-1}}{h_J}\left( \frac{a(x_{J-1/2})h_{J}}{2h_{-}} - \frac{d(x_{J-1/2})}{h_{-}} \right) + \frac{w_{J}}{h_J} \left( \frac{a(x_{J-1/2}) h_{J-1}}{2h_{-}} + \frac{d(x_{J-1/2})}{h_{-}} \right) - \frac{g_{R}(x_R)}{h_J}

The :math:`\theta`-method
*************************

The :math:`\theta`-method is an approach which improves the stability and numerical accuracy when integrating a partial differential equation in time. It consists of writing the equation as the time average of the current and future time step. When :math:`\theta=0` a fully explicit scheme is recovered in which the future state of the system is derived purely from the current state. Conversely,  :math:`\theta=1` gives a fully implicit formalism, in which a linear system of equation is solved to determine the future state. Setting :math:`\theta=1/2` results in an average of these two limits and it is generally causes the the Crank-Nicolson method. Crank-Nicolson provides unconditionally stable iterations for the advection and diffusion equations, and the improves the time integration (it corresponds to a trapezium integration in the time domain).


.. math::
	\frac{w_j^{n+1} - w_j^n}{k} = \theta r_a w_{j-1}^{n+1} + (1 - \theta) r_a w_{j-1}^{n} + \theta r_b w_{j-1}^{n+1}  + (1 - \theta) r_b w_{j-1}^{n} +  \theta r_c w_{j-1}^{n+1}  + (1 - \theta) r_c w_{j-1}^{n} + \bar{s}_j^n

In the above :math:`k` stands for the difference in time, the :math:`n+1` are the terms at the future time point, and :math:`n` terms are the current time point. The ":math:`r`" terms are the coefficients of the semi-discretised equation. Moving the unknowns to the left hand side,

.. math::
	w_j^{n+1} - \theta k r_a w_{j-1}^{n+1} - \theta k r_b w_{j}^{n+1} - \theta k r_c w_{j+1}^{n+1} = w_j^{n} + (1 - \theta) k r_a w_{j-1}^{n} + (1 - \theta) k r_b w_{j}^{n} + (1 - \theta) k r_c w_{j1}^{n} + k \bar{s}_j^n
	
Defining the coefficients for the interior,

.. math::
	r_a & = \frac{k}{h_j} \left( \frac{a(x_{j-\frac{1}{2}}) h_j}{2h_{-}}  - \frac{d(x_{j-\frac{1}{2}})}{h_{-}}\right) \\
	r_b & = \frac{k}{h_j}\left( \frac{a(x_{j-\frac{1}{2}})h_{j-1}}{2h_{-}} - \frac{a(x_{j+\frac{1}{2}})h_{j+1}}{2h_{+}} + 	\frac{d(x_{j-\frac{1}{2}})}{h_{-}} + \frac{d(x_{j+\frac{1}{2}})}{h_{+}}  \right)\\
	r_c & = \frac{k}{h_j} \left( \frac{-a(x_{j+\frac{1}{2}})h_j}{2h_{+}} - \frac{d(x_{j+\frac{1}{2}})}{h_{+}} \right)

the left hand side boundary,

.. math::
	\alpha_b = \frac{k}{h_1}\left( -\frac{a(x_{3/2})h_2}{2h_{+}} + \frac{d(x_{3/2})}{h_{+}} \right) \\
	\alpha_c = \frac{k}{h_1}\left( -\frac{a(x_{3/2})h_1}{2h_{+}} - \frac{d(x_{3/2})}{h_{+}} \right)

and the right hand side boundary,

.. math::
	\beta_a = \frac{k}{h_J}\left( \frac{a(x_{J-1/2})h_J}{2h_{-}} - \frac{d(x_{J-1/2})}{h_{-}} \right) \\
	\beta_b = \frac{k}{h_J}\left( \frac{a(x_{J-1/2})h_{J-1}}{2h_{-}} + \frac{d(x_{J-1/2})}{h_{-}} \right)


The linear system including becomes,

.. math::
	\begin{align} 
	\begin{pmatrix}
	1-\theta \alpha_b & -\theta \alpha_c    &        &       & 0   \\
	-\theta r_a & 1-\theta r_b    & -\theta r_c    &       &     \\
	    & \ddots & \ddots & \ddots&     \\
	    &        &  - \theta r_a   & 1-\theta r_b   & -\theta r_c \\
	 0  &        &        & -\theta \beta_a   & 1-\theta \beta_b
	\end{pmatrix}
	\begin{pmatrix}
	    w_1^{n+1} \\
	    w_2^{n+1} \\
	    \vdots \\
	    w_{J-1}^{n+1} \\
	    w_J^{n+1} \\
	\end{pmatrix} = 
	\begin{pmatrix}
	1+(1-\theta)\alpha_b & (1-\theta)\alpha_c    &        &       & 0   \\
	(1-\theta)r_a & 1+(1-\theta)r_b    & (1-\theta)r_c    &       &     \\
	    & \ddots & \ddots & \ddots&     \\
	    &        &  (1-\theta)r_a   & 1+(1-\theta)r_b   & (1-\theta)r_c \\
	 0  &        &        & (1-\theta)\beta_a   & 1+(1-\theta)\beta_b
	\end{pmatrix}
	\begin{pmatrix}
	    w_1^n \\
	    w_2^n \\
	    \vdots \\
	    w_{J-1}^n \\
	    w_J^n \\
	\end{pmatrix} + k
	\begin{pmatrix}
	    s_1^n + g_R(x_L)/h_1 \\
	    s_2^n \\
	    \vdots \\
	    s_{J-1}^n \\
	    s_J^n - g_R(x_R)/h_J\\
	\end{pmatrix}
	\end{align}

Notice how the boundary conditions appear in the vector on the right hand side.



Aside :math:`-` Linear interpolation between cell centre and face values
=========================================================================

In general, linear interpolation between two points :math:`(x_0, x_1)` can be used to find the value of a function at :math:`f(x)`,

.. math::
	f(x) = \frac{x - x_1}{x_0 - x_1}f(x_0) + \frac{x - x_0}{x_1 - x_0}f(x_1)

In a cell centred grid we know the value of the variable :math:`w` at difference points, :math:`w_j` and :math:`w_{j+1}`. We can apply the linear interpolation formulae above to determine value at cell face :math:`w_{j+1/2}`.

.. math::
	w_{j+1/2} =  \frac{x_{j+1/2} - x_{j+1}}{x_{j} - x_{j+1}} w_j + \frac{x_{j+1/2} - x_j}{x_{j+1} - x_j} w_{j+1} 

This can be simplified firstly by using function to represent the distance between cell centres,

.. math::
	h_{-} = x_j - x_{j-1} \quad h_{+} = x_{j+1} - x_{j}

to give, 

.. math::
	w_{j+1/2} = \frac{x_{j+1} - x_{j+1/2}}{h_{+}} w_j + \frac{x_{j+1/2} - x_j}{h_{+}} w_{j+1}

This expression still contains :math:`x_{j+1/2}` which we can simplify further by using an expression for the position of cell centres,

.. math::
	x_j = \frac{1}{2} \left( x_{j-\frac{1}{2}} + x_{j+\frac{1}{2}} \right) \quad x_{j+1} = \frac{1}{2} \left( x_{j+\frac{1}{2}} + x_{j+\frac{3}{2}} \right)


Note, this expression is still valid of non-uniform grids, it simply says that cell centres are always equidistant from two faces. Rearranging the above expression and substituting in for :math:`x_{j}` and :math:`x_{j+1}` terms gives, 

.. math::
	w_{j+1/2} = \frac{\frac{1}{2} \left( x_{j+\frac{1}{2}} + x_{j+\frac{3}{2}} \right) - x_{j+1/2}}{h_{+}} w_j + \frac{x_{j+1/2} - \frac{1}{2} \left( 	x_{j-\frac{1}{2}} + x_{j+\frac{1}{2}} \right)}{h_{+}} w_{j+1}


Finally, by defining the distance between vertices as, :math:`h_j = x_{j+\frac{1}{2}} - x_{j-\frac{1}{2}}`, we can simplify to the following expression,

.. math::
	w_{j+1/2} = \frac{h_{j+1}}{2h_{+}} w_j + \frac{h_j}{2h_{+}} w_{j+1}


Similarly the :math:`w_{j-1/2}` can be found,

.. math::
	w_{j-1/2} = \frac{h_{j}}{2h_{-}} w_{j-1} + \frac{h_{j-1}}{2h_{-}} w_{j}

