---
layout: post
title:  "Fractional-Step Pressure Projection Algorithm on Staggered Grid for Incompressible Navier-Stokes"
date:   2024-03-28 15:00:00 +0700
categories: cfd
---

# Incompressible Navier-Stokes Equations with Newtonian Fluid and Constant Viscosity

In Gibbs notation:

$$
\nabla \cdot \textbf{u} = 0
$$

$$\frac{\partial \textbf{u}}{\partial t} = -\nabla \cdot \{\textbf{u}\otimes\textbf{u}\} - \nabla P + \nu \nabla^2 \textbf{u}$$

In tensor notation:

$$
\frac{\partial u_i}{\partial x_i} = 0
$$

$$ \frac{\partial u_i}{\partial t} = -\frac{\partial(u_i u_j)}{\partial x_j} - \frac{\partial P}{\partial x_i} + \nu \frac{\partial^2 u_i}{\partial x_j \partial x_j}$$

Note that $$P = p/\rho$$ and is called the *kinematic pressure*.

# Fractional-Step Pressure Projection Formulation

At any timestep $$n$$, we march forward in time to get the velocity and pressure field at the next timestep $$n+1$$. First the intermediate velocity field $$\textbf{u}^F$$ is obtained from the current velocity field $$\textbf{u}^n$$:

$$
\begin{equation}
\textbf{u}^F = \textbf{u}^n + \Delta t (-\nabla \cdot \{\textbf{u}\otimes\textbf{u}\} + \nu \nabla^2 \textbf{u})^n
\end{equation}
$$

Notice that this intermediate velocity field is calculated without taking into account the pressure gradient. Then $$\textbf{u}^F$$ is corrected to get $$\textbf{u}^{n+1}$$:

$$
\begin{equation}
\textbf{u}^{n+1} = \textbf{u}^F - \Delta t \nabla P^{n+1}
\end{equation}
$$

where $$\phi$$ is obtained by solving the Poisson equation:

$$
\begin{equation}
\nabla^2 P^{n+1} = \frac{1}{\Delta t} \nabla \cdot \textbf{u}^F
\end{equation}
$$

This equation is the result of taking the divergence of the velocity field correction equation and enforcing continuity on the velocity field at the next timestep ($$\nabla \cdot \textbf{u}^{n+1} = 0$$). The pressure $$P$$ is seen as a scalar field that forces $$\textbf{u}^{n+1}$$ to be solenoidal. This is the role pressure plays in numerical incompressible Navier-Stokes when the pressure correction formulation is being used.

Note also that the above equations are written under the assumption that we are marching forward in time in an explicit Euler manner. The actual equations depend on how time marching is done. For simplicity, let us stick to forward Euler for now. We will discretize the equations with finite-differencing on an internal cell first then deal with boundary conditions later.

# Discretization of the Advection Term
The advection term in Equation 1 will now be discretized. We first consider the first velocity component which for brevity will now be called just $$u$$:

$$
\frac{\partial(u_1 u_j)}{\partial x_j} = \frac{\partial(u_1 u_1)}{\partial x_1} + \frac{\partial(u_1 u_2)}{\partial x_2} = \frac{\partial(u^2)}{\partial x} + \frac{\partial(uv)}{\partial y}
$$

Our job now is to express these partial derivatives with the velocity information in or staggered grid. At an internal grid point $$(i,j)$$ where $$u$$ is stored, the central-difference expressions will now be written. Refer to FIGURE for the compass notation ($$\text{N-E-S-W}$$) and where things are.

$$
\frac{\partial(u^2)}{\partial x} \approx \frac{u_\text{E}^2 - u_\text{W}^2}{x_{i+1} - x_i}
$$

The values at $$\text{E}$$ and $$\text{W}$$ are approximated with linear interpolation:

$$
\begin{aligned}
u_\text{E} &\approx (1-\xi_\text{E})u_{i,j} + \xi_\text{E}u_{i+1,j} \\
u_\text{W} &\approx (1-\xi_\text{W})u_{i-1,j} + \xi_\text{W}u_{i,j}
\end{aligned}
$$

and the interpolation fractions are

$$
\begin{aligned}
\xi_\text{E} &= \frac{x_{i+1} - x_i}{x_{i+2}-x_i}\\ 
\xi_\text{W} &= \frac{x_i - x_{i-1}}{x_{i+1}-x_{i-1}}\\ 
\end{aligned}
$$

These fractions measure how close to point being interpolated for is to the first (leftmost or downmost) point. The value of zero means coincidence to the first point and the value of one to the second.

Now for the next term:

$$
\frac{\partial(uv)}{\partial y} \approx \frac{u_\text{N} v_\text{N} - u_\text{S} v_\text{S}}{\frac{1}{2}(y_{j+1} - y_{j-1})}
$$

The values at $$\text{N}$$ and $$\text{S}$$ are approximated with linear interpolation:

$$
\begin{aligned}
u_\text{N} &\approx \frac{1}{2}(u_{i,j} + u_{i,j+1}) \\
u_\text{S} &\approx \frac{1}{2}(u_{i,j-1} + u_{i,j}) \\
v_\text{N} &\approx \frac{1}{2}(v_{i,j} + v_{i+1,j}) \\
v_\text{S} &\approx \frac{1}{2}(v_{i,j-1} + v_{i+1,j-1})
\end{aligned}
$$

Similarly for the second velocity component:

$$
\frac{\partial(u_2 u_j)}{\partial x_j} = \frac{\partial(u_2 u_1)}{\partial x_1} + \frac{\partial(u_2 u_2)}{\partial x_2} = \frac{\partial(uv)}{\partial x} + \frac{\partial(v^2)}{\partial y}
$$