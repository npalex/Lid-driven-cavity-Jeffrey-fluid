# **Lid-driven cavity flow of an incompressible Jeffrey fluid**

&emsp; This program solves the continuity, incompressible Cauchy momentum equations, and Jeffrey's viscoelastic fluid model in 2D, given by

$$ \frac{\partial u}{\partial x} + \frac{\partial v}{\partial y} = 0, $$

$$ Re\left( \frac{\partial u}{\partial t} + \frac{\partial u^2}{\partial x} + \frac{\partial (u v)}{\partial y} 
	+ \frac{\partial p}{\partial x}\right) = \frac{\partial \tau_{xx}}{\partial x} + \frac{\partial \tau_{xy}}{\partial y} 
    + \beta \left( \frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} \right), $$

$$ Re\left( \frac{\partial v}{\partial t} + \frac{\partial (u v)}{\partial x} + \frac{\partial v^2}{\partial y} 
	+ \frac{\partial p}{\partial y}\right) = \frac{\partial \tau_{xy}}{\partial x} + \frac{\partial \tau_{yy}}{\partial y} 
    + \beta \left( \frac{\partial^2 v}{\partial x^2} + \frac{\partial^2 v}{\partial y^2} \right), $$
	
$$ Wi\frac{\partial \tau_{xx}}{\partial t} + \tau_{xx} = 2\frac{\partial u}{\partial x}, $$

$$  Wi\frac{\partial \tau_{xy}}{\partial t} + \tau_{xy} = \frac{\partial v}{\partial x} + \frac{\partial u}{\partial y}, $$

and

$$ Wi\frac{\partial \tau_{yy}}{\partial t} + \tau_{yy} = 2\frac{\partial v}{\partial y}. $$

where $u$ and $v$ are the components of the local fluid velocity in the $x$ and $y$ directions, $\tau_{xy}$ is the local "polymeric" shear stress, 
and $\tau_{xx}$ and $\tau_{yy}$ are local, polymeric normal stresses.

&emsp; The Reynolds number $Re$, Weissenberg number $Wi$, and viscosity ratio $\beta$ are defined according to 

$$Re = \frac{\rho U L}{\eta},$$

$$Wi = \lambda \frac{U}{L},$$

and

$$ \beta = \frac{\eta_s}{\eta}. $$

Here, $\rho$, $\eta$, $\eta_s$, $U$, and $L$ are the over all fluid density, fluid viscosity, solvent viscosity, lid speed, and lid length, respectively. The time scale $\bar t$ 
for this problem is defined by the characteristic shear rate according to $\bar t = \frac{L}{U}$ and the inertial pressure scale was chosen, 
equal to $\bar p = \rho U^2$. The boundary conditions are no-slip and no-flow at the cavity walls and the fluid is initially at rest. 
 
## **Numerical Scheme:**
&emsp; Following Lee and Leveque (2003),<sup>1</sup> a fractional step approach is used to solve the system of equations above, which decomposes the problem into the following steps:


### **Step 1. Solve the convection equation:**

&emsp;A Roe solver $\textemdash$ i.e., a locally linear, approximate Riemann solver based on the Godunov method $\textemdash$ 
is employed to evaulate $q = (u,v,\tau_{xx},\tau_{xy},\tau_{yy})$ satisfying:

$$ q_t + \hat A \cdot q_x + \hat B \cdot q_y = 0 . $$

where the Roe matrices $\hat A$ and $\hat B$ are approximate Jacobian matrices given by

$$ \hat A =        \begin{bmatrix} 
                                2\hat u & 0 & -\frac{1}{Re} & 0 & 0 \\
								\hat v & \hat u & 0 & -\frac{1}{Re} & 0 \\
								-\frac{2}{Wi} & 0 & 0 & 0 & 0 \\
								0 & -\frac{1}{Wi} & 0 & 0 & 0 \\
                                0 & 0 & 0 & 0 & 0 \end{bmatrix}, $$
				
and

$$ \hat B =        \begin{bmatrix} 
                                \hat v & \hat u & 0 & -\frac{1}{Re} & 0 \\
								0 & 2\hat v & 0 & 0 & -\frac{1}{Re} \\
								0 & 0 & 0 & 0 & 0 \\
								-\frac{1}{Wi} & 0 & 0 & 0 & 0 \\
                                0 & -\frac{2}{Wi} & 0 & 0 & 0 \end{bmatrix}, $$

Here, $\hat u$ and $\hat v$ are Roe averages (in this case, they are linear interpolations of cell-centered velocities) defined at the edge of each grid cell. For example, the Roe averages used to evaluate the matrix $A$ are given by

$$ \hat u_{i-\frac{1}{2},j} = \frac{U_{i,j} + U_{i-1,j}}{2}$$

and

$$ \hat v_{i-\frac{1}{2},j} = \frac{V_{i,j} + V_{i-1,j}}{2}$$

Dimensional splitting via the donor cell upwind method (DCU) is used to advanced the cell-centered velocities and stresses $Q=(U,V,T_{xx},T_{xy},T_{yy})$ forward in time via sweeps in the x-direction

$$Q_{i,j}^{\*} = Q^n_{i,j} - \frac{\Delta t}{\Delta x} \left( F_{i+\frac{1}{2},j}^{n} - F_{i-\frac{1}{2},j}^{n}\right). $$

followed by sweeps in the y-direction

$$Q_{i,j}^{\*\*} = Q_{i,j}^{\*} - \frac{\Delta t}{\Delta y} \left( G_{i,j+\frac{1}{2}}^{\*} - G_{i,j-\frac{1}{2}}^{\*}\right). $$

where $F_{i-\frac{1}{2},j}$ is the numerical flux at the interface between cells $(i,j)$ and $(i-1,j)$ for the 1-dimensional problem in the x-direction and, similarly, $G_{i,j-\frac{1}{2}}$ is the flux at the interface between cells $(i,j)$ and $(i,j-1)$ for the 1D problem in the y-direction. In addition, monotenzied central flux limiters are used to achieve second order accuracy for this step where the solution is smooth. 

### **Step 2. Solve the diffusion equation:** 

During this step, the following equation is solved

$$ q_t = \psi(q), $$

where the source terms are given by

$$ \psi(q) = \begin{bmatrix} 
                                \frac{\beta}{Re}(u_{xx} + u_{yy}) \\
								\frac{\beta}{Re}(v_{xx} + v_{yy}) \\
								-\frac{1}{Wi}\tau_{xx} \\
								-\frac{1}{Wi}\tau_{xy} \\
                                -\frac{1}{Wi}\tau_{yy} \end{bmatrix}. $$
								
&emsp; An alternating direction implicit (ADI) method is employed to update the velocities for diffusion. 
Two difference equations are used to advance successive time steps of $\frac{\Delta t}{2}$. The first equation for the u-velocity is implict in the x-direction

$$ U_{i,j}^{\*\*\*} = U_{i,j}^{\*\*} + \frac{\alpha}{2}\left(\delta^2_x U^{\*\*\*} + \delta^2_y U^{\*\*}\right) $$

and the second equation for the u-velocity is implicit in the y-direction

$$ \widetilde U_{i,j} = U_{i,j}^{\*\*\*} + \frac{\alpha}{2}\left(\delta^2_x U^{\*\*\*} + \delta^2_y \widetilde U\right). $$

The parameter $\alpha$ is defined by

$$ \alpha = \frac{\beta}{Re}\frac{\Delta t}{(\Delta x)^2},$$

and $\delta^2_x$ denotes the central difference of the 2nd partial derivative with respect to $x$. Similar equations are used to update the v-velocity. 

### **Step 3. Update the edge velocities for diffusion:**

&emsp; The velocities at the edges of each grid cell are updated via linear interpolation according to

$$ \widetilde u_{i-\frac{1}{2},j} = \frac{\widetilde U_{i-1,j} + \widetilde U_{i,j}}{2} $$

and

$$ \widetilde v_{i-\frac{1}{2},j} = \frac{\widetilde V_{i-1,j} + \widetilde V_{i,j}}{2} $$

### **Step 4. Compute the pressure distribution:**

&emsp;So far, the velocity field $(\widetilde u, \widetilde v)$ is not divergence free. In order to satisfy continuity, $(\widetilde u, \widetilde v)$ is projected into a divergence-free vector field by correcting the result for pressure-driven flow via

$$ u_{i,j}^{n+1} = \widetilde u_{i,j} -\Delta t\nabla p^{n+1}$$

and 

$$ v_{i,j}^{n+1} = \widetilde v_{i,j} -\Delta t\nabla p^{n+1}$$

&emsp; In order to acquire the pressure distribution, the divergence of the equation above provides a Laplacian equation for the pressure,

$$ \nabla^2 p^{n+1} = \frac{1}{\Delta t} \left( \frac{\partial \widetilde u_{i,j}}{\partial x} + \frac{\partial \widetilde v_{i,j}}{\partial y} \right) ,$$

which is then discretized with Neumann boundary conditions to produce a system of linear equations, $Ax = b$. 
However, the matrix $A$ is singular because the equation set has an inifinite number of solutions within an arbitrary reference pressure. 
Hence, a ficticious source term $C_0 p^{n+1}$ has been added with proportionality constant $C_0$, 
which is defined on the order of 1e-9 to render the influence of the source negligble, so that $A$ is non-singular.  

### **Step 5. Update the edge velocities for pressure-driven flow:**

&emsp; The edge velocities are determined by using central differences for the pressure gradient via

$$ u_{i-\frac{1}{2},j}^{n+1} = \widetilde u_{i-\frac{1}{2},j} - \Delta t \left( \frac{p_{i,j}^{n+1} - p_{i-1,j}^{n+1}}{\Delta x}\right) $$

and

$$ v_{i,j-\frac{1}{2}}^{n+1} = \widetilde v_{i,j-\frac{1}{2}} - \Delta t \left( \frac{p_{i,j}^{n+1} - p_{i,j-1}^{n+1}}{\Delta y}\right) $$

### **Step 6. Update the cell-centered velocities for pressure-driven flow:**

&emsp;Finally, the cell-centered velocities are determined by using central differences for the pressure gradient according to

$$ U_{i,j}^{n+1} = \widetilde U_{i,j} - \Delta t \left( \frac{p_{i+1,j}^{n+1} - p_{i-1,j}^{n+1}}{2 \Delta x}\right) $$

and

$$ V_{i,j}^{n+1} = \widetilde V_{i,j} - \Delta t \left( \frac{p_{i,j+1}^{n+1} - p_{i,j-1}^{n+1})}{2 \Delta y}\right), $$


## **Results**:

### **Velocity and pressure distribution for $Re = 0.1$, $Wi = 0.1$, and $\beta = 0$ with $CFL \leq 0.5$ using the DCU method on a 51x51 cell grid over the time interval $[0, 1]$**

### **Streamlines and normal stress distribution for $Re = 0.1$, $Wi = 0.1$, and $\beta = 0$ with $CFL \leq 0.5$ using the DCU method on a 51x51 cell grid over the time interval $[0, 1]$**

## **References**:

1.	L.Lee and R.J.LeVeque, 2003. An immersed interface method for incompressible Navier
		Stokes equations. SIAM J. Sci. Comput., 25, 832–856.

2.	Clawpack Development Team (2023), Clawpack Version 5.9.2,
		http://www.clawpack.org, doi: 10.5281/zenodo.10076317

3.	R. J. LeVeque, 1997. Wave propagation algorithms for multi-dimensional 
		hyperbolic systems. J. Comput. Phys. 131, 327–353.

4.	R. J. LeVeque. Finite Volume Methods for Hyperbolic Problems. Cambridge 
		University Press, Cambridge, UK, 2002.