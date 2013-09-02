Solving linear equations
------------------------

Provided the equation is linear, meaning that the coefficients, nor the reaction term depend on the solution variable a time-stepping approach can be used to solve the equation.

Time stepping the advection-diffusion equation
**********************************************

Here we will rearrange the :math:`\theta`-method form of the advection-diffusion equation,

.. math::
	\frac{w^{n+1} - w^{n}}{\tau} & =  \theta F(w^{n+1}) + (1-\theta)F(w^{n})
	
into the form of a linear system :math:`A\cdot x = d`. Firstly lets move all terms involving future time points the l.h.s,

.. math::
    w^{n+1} - \theta \tau F(w^{n+1}) = w^{n} + (1-\theta) F(w^{n+1})

Replacing the :math:`F` terms with the full matrix expressions and factoring yields,

.. math::
    \underbrace{(I - \theta\tau M^{n+1})}_{A}\underbrace{w^{n+1}}_{x} = \underbrace{(I + (1-\theta)\tau M^{n})w^{n}}_{d}

Time-stepping involves solving this equation iteratively. First the initial conditions :math:`w^0` is used to calculate the solution variable at the next point in time :math:`w^1`, then values of the solution variable are updated such that :math:`w^1` is used to calculate :math:`w^2`, and so on and so forth. This procedure is shown in the Figure below,

.. figure:: img/time_stepping_linear_advection_diffusion.png
   :scale: 50 %
   :alt: Time stepping the linear advection-diffusion equation.
   :align: center

.. note::
	Parameters used for this simulation. :math:`a=1`, :math:`d=1\times 10^{-3}`, :math:`h=5 \times 10^{-3}`, :math:`\tau=5 \times 10^{-4}` with exponentially fitted discretisation. The initial and boundary conditions where, :math:`u(x,t_{0}) =sin(\pi x)^{100}` with :math:`u(x_1)=1` and :math:`u_x(x_J)=0`. The equation was time stepped for 400 iterations, the plots so the solution at times :math:`t=0,0.05,0.1,0.15,0.2`.