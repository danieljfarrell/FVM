Example calculations
--------------------

.. note::
	For all simulations :math:`a=1`, :math:`d=0.001` on a domain where :math:`0\leq x \leq 1`, with mesh spacing :math:`h=0.02` and time step :math:`\tau=0.01`. The equations are solver with the boundary conditions :math:`u(0,t)=1` and :math:`u(1,t)=0` and initial conditions :math:`u(x,0)=sin(\pi x)^{100}`. With these boundary conditions a boundary layer will form near :math:`x=1`.

Uniform grid
************

First we test the finite-volume method using a standard uniform grid. Note that the Peclet number for the above parameters is :math:`\mu=20` so the central discretisation scheme is not stable as illustrated by the oscillations in the solution. The upwind scheme does not have a stability criteria related to the Peclet number so the solutions for the *upwind* case are smooth. Finally, the *exponentially fitted* scheme has automatically weighted in favour of the upwind discretisation, the value of :math:`\kappa\approx 0.9`.

.. raw:: html

    <div style="margin-top:10px;">
	  <iframe src="http://player.vimeo.com/video/69527955" width="480" height="460" frameborder="0" webkitAllowFullScreen mozallowfullscreen allowFullScreen></iframe>
    </div>

Random grid
***********

Although a random grid is of no practical use it is a good test of the code because bugs are more likely to show up when symmetry has been reduced. I the follow simulations :math:`\text{min}(\kappa)` =0.004 and :math:`\text{max}(\kappa)` =0.98 so the discretisation scheme is abruptly changing from cell to cell.

.. raw:: html

    <div style="margin-top:10px;">
	  <iframe src="http://player.vimeo.com/video/69528243" width="480" height="460" frameborder="0" webkitAllowFullScreen mozallowfullscreen allowFullScreen></iframe>
    </div>

Nonuniform grid
***************

Nonuniform grids can be used to reduce to *improve* the solution as shown here. The following simulation contains the same number of cells as the previous simulations however the cell centres are clustered towards the right hand boundary. The increased density of cells allow the boundary layer to be resolved. The transient solution computed from the *central* scheme is still significantly affected by the high Peclet number but it is interesting to observe that the steady-state solution of all three methods are very similar. Furthermore, :math:`\kappa\approx0.01-0.1` in the region of the boundary layer which implies the local value of Peclet number as been reduced enough so allow the exponentially fitted scheme to weight in favour of the higher accuracy central discretisation. 

.. raw:: html

    <div style="margin-top:10px;">
	  <iframe src="http://player.vimeo.com/video/69528242" width="480" height="460" frameborder="0" webkitAllowFullScreen mozallowfullscreen allowFullScreen></iframe>
    </div>
