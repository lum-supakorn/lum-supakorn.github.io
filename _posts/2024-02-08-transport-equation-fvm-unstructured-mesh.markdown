---
layout: post
title:  "Notes on Discretization of Scalar Transport Equation with Finite-Volume Method on Unstructured Mesh"
date:   2024-02-08 14:04:56 +0700
categories: cfd
---
## Why?
I wrote this note to keep track of all the equations while researching finite-volume methods. I have derived all these equations in a notebook and found myself having to re-derive them for a couple of times now because I just could not find where I wrote what I needed!

I am trying to solve a transport problem on an unstructured mesh and slowly build this up to a [RANS][rans] solver. This problem is from section 4.7.2 of [Computational Methods for Fluid Dynamics][ferziger] except that in the book it was a structured mesh. The discretization method is also from that book with a bit of help from [The Finite Volume Method in Computational Fluid Dynamics][moukalled]. You can find much more information on finite-volume discretization there. This post is intended to be more of my personal note and less of an educational material.

Please note that I have not made any meaningful effort in analyzing the errors introduced by a series of approximation methods shown in this note. I am aiming for second-order accuracy and shall present the analysis soon.

In the end I want a linear system:

$$ \textbf{A}\bm{\phi} = \textbf{b} $$

where $$\bm{\phi}$$ is a vector of unknown scalars, each defined at the position of each cell in the field. After getting this system I can analyze it and see how it can be solved with parallel, throughput-intensive algorithms. This note describes the process of getting $$\textbf{A}$$ and $$\textbf{b}$$.

## Scalar Transport Equation

The equation we are working with is the steady-state scalar transport (convection-diffusion) equation with no sources. The equation is integrated over an arbitrary volume (cell), and the divergence theorem is applied to yield surface integrals:

$$ \int_S\rho\phi\textbf{v}\cdot \textbf{n}\mathrm{d}S = \int_S \Gamma \nabla \phi\cdot\textbf{n}\mathrm{d}S $$

Approximating the surface integral with the midpoint rule on $$k$$ cell faces yields:

$$ \sum_k \rho \phi_k \textbf{v}_k \cdot \textbf{S}_k = \sum_k \Gamma (\nabla \phi)_k \cdot \textbf{S}_k $$

where the subscript $$k$$ indicates the value at the midpoint of each face.

To build the linear system, we figure out the elements of matrix $$\textbf{A}$$ and vector $$\textbf{b}$$. Row $$i$$ of $$\textbf{A}$$ and element $$i$$ of $$\textbf{b}$$ represents the discretized scalar transport equation when it is applied on cell $$i$$ of the mesh. The equation for cell $$i$$ is written as a linear combination of the unknown scalar values $$\phi$$ at the center of that particular cell $$\mathrm{P}$$ its the neighbor cells $$\mathrm{N}_k$$. Element $$\textbf{A}(i,j)$$ is the coefficient associated with the unknown of cell $$j$$ in the equation for cell $$i$$. $$\textbf{b}(i)$$ collects the constant terms and is the right-hand side of equation $$i$$.

All of the equations in the resulting linear system would look like this:

$$
\left(a_{\mathrm{P}}\phi_{\mathrm{P}} + \sum_k a_{\mathrm{N}_k}\phi_{\mathrm{N}_k}\right)_i = \mathrm{RHS}_i
$$

and it is expanded from:

$$
F^{\mathrm{c}} - F^{\mathrm{d}} = 0
$$

where $$F^{\mathrm{c}}$$ is the total convection flux and $$F^{\mathrm{d}}$$ is the total diffusion flux. I must be careful to remember the minus sign in front of the diffusion term when collecting the coefficients.

Building $$\textbf{A}$$ and $$\textbf{b}$$ is just a matter of calculating for these coefficients and the right-hand side. This depends on the geometry of and the boundary conditions imposed on each cell. Finally, a coefficient is a combination of one from convection and the other from diffusion:

$$
a = a^\mathrm{c} + a^\mathrm{d}
$$

## Convection Term
We can now write the discretized equation in terms of the unknown $$\phi$$'s which is defined at the center of each cell. Let's start with the convection summation first:

$$
F^{\mathrm{c}} = \sum_k \rho \phi_k \textbf{v}_k \cdot \textbf{S}_k
$$

We denote each term as:

$$ F^{\mathrm{c}}_k = \rho \phi_k \textbf{v}_k \cdot \textbf{S}_k $$

Since the velocity field and the flow density are known and mesh information is available, $$\rho$$, $$\textbf{v}_k$$, and $$\textbf{S}_k$$ are known. The only thing left to do is to express $$\phi_k$$, the scalar value at the midpoint of each face, in terms of $$\phi$$ at the center of each cells, our unknowns in the upcoming linear system. This depends on how we express the face-centered value in terms of the cell-centered value and the boundary condition imposed on each face.

# Convection Term: Internal Faces
A linear interpolation is applied to express $$\phi_k$$ of an internal face as a linear combination of the scalar of the two straddling cells, the owner cell $$\mathrm{P}$$ and the neighbor cell $$\mathrm{N}_k$$. Let's call the result of this first approximation $$\phi_{k'}$$:

$$
\phi_{k'} = (1-\xi_k)\phi_{\mathrm{P}} + \xi_k\phi_{\mathrm{N}_k}
$$

where $$ 0 \leq \xi_k \leq 1$$ is the interpolation factor. It tells us how close the face midpoint is to the two cell centers. A value of 0 means that the face center is at the center of cell $$\mathrm{P}$$ while a value of 1 means that the face center is at the center of cell $$\mathrm{N}_k$$.

This interpolation factor is calculated by projecting the displacement vector $$\bm{r}_k - \bm{r}_{\mathrm{P}}$$ onto $$\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}$$ and get the relative magnitude:

$$
\xi_k = \frac{(\bm{r}_k - \bm{r}_{\mathrm{P}}) \cdot (\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}})}{|\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}|^2}
$$

But that's not all. The face midpoint where $$\phi_k$$ lies might not be on the segment connecting the centers of cell $$\mathrm{P}$$ and $$\mathrm{N}_k$$, so another approximation is needed. The figure below illustrates this:

|![Internal face of finite-volume cell](/images/transport-equation-fvm-unstructured-mesh/fvm_cell_internal_faces.svg)
|:--:| 
| *Figure 1:* The center of a face $$k$$ may not lie on the segment connecting the cell centers |

We address this by introducing an auxiliary point $$k'$$. In fact, the $$\phi_{k'}$$ introduced earlier is the $$\phi_k$$ defined on this auxiliary point. It is the point which lies on the segment joining the straddling cell centers. $$\textbf{r}_{k'}$$ is obtained by projecting $$\bm{r}_k - \bm{r}_{\mathrm{P}}$$ and $$\bm{r}_{\mathrm{N}_k} - \bm{r}_k$$ onto $$\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}$$. Whichever projected vector is shorter gets to point to $$k'$$. The figure below illustrates the case where $$\bm{r}_k - \bm{r}_{\mathrm{P}}$$ is shorter and thus selected.

|![Internal face of finite-volume cell correction](/images/transport-equation-fvm-unstructured-mesh/fvm_cell_internal_faces_correction.svg)
|:--:| 
| *Figure 2:* The auxiliary point $$k'$$ |

The position of $$k'$$ is therefore:

$$
\textbf{r}_{k'} = \textbf{r}_{\mathrm{P}} + \alpha(\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}})
$$

where

$$
\alpha = \mathrm{min}\left(\frac{(\bm{r}_k - \bm{r}_{\mathrm{P}})\cdot(\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}})}{|\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}|^2}, \frac{(\bm{r}_{\mathrm{N}_k} - \bm{r}_k)\cdot(\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}})}{|\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}|^2}\right)
$$

The two approximations are combined into the following expression for $$\phi_k$$:

$$
\phi_k = \phi_{k'} + (\nabla \phi)_{k'} \cdot (\bm{r}_k-\bm{r}_{k'})
$$

The second term describes the change in $$\phi_k$$ as we move from $$\bm{r}_{k'}$$ to $$\bm{r}_k$$ If information on $$\nabla \phi$$ is available at the center of each cell, $$(\nabla \phi)_{k'}$$ can be obtained with linear interpolation in the same manner as $$\phi_{k'}$$. Since $$\nabla \phi$$ depends on $$\phi$$ on each cell and vice versa, an iterative method is needed to obtain both. 

One possible way is to first omit this correction altogether on the first iteration, construct the linear system without this correction, solve for $$\bm{\phi}$$, determine $$\nabla \phi$$, apply the correction, and repeat until convergence. At iteration $$i$$, the correction calculated from the scalars of the previous iteration $$i-1$$ is used. This is called *deferred correction*.

<div class="post-blue">
The convection flux for internal faces is thus:

$$
F^{\mathrm{c}}_k = \rho \{\overbrace{(1-\xi_k)\phi_{\mathrm{P}} + \xi_k\phi_{\mathrm{N}_k}}^{\phi_{k'}} + \overbrace{(\nabla \phi)_{k'} \cdot (\bm{r}_k-\bm{r}_{k'})}^\text{correction}\} (\textbf{v}_k \cdot \textbf{S}_k)
$$

Collecting the coefficients:
$$
F^{\mathrm{c}}_k = \{\rho(1-\xi_k)(\textbf{v}_k \cdot \textbf{S}_k)\}\phi_{\mathrm{P}} + \{\rho\xi_k(\textbf{v}_k \cdot \textbf{S}_k)\}\phi_{\mathrm{N}_k} + \rho(\textbf{v}_k \cdot \textbf{S}_k)(\nabla \phi)_{k'} \cdot (\bm{r}_k-\bm{r}_{k'})
$$

Its contribution to the linear system:

$$
\textbf{A}(i,i) \mathrel{+}= \rho(1-\xi_k)(\textbf{v}_k \cdot \textbf{S}_k)
$$

$$
\textbf{A}(i,j) \mathrel{+}= \rho\xi_k(\textbf{v}_k \cdot \textbf{S}_k)
$$

$$
\textbf{b}(i) \mathrel{+}= -\rho(\textbf{v}_k \cdot \textbf{S}_k)(\nabla \phi)_{k'} \cdot (\bm{r}_k-\bm{r}_{k'})
$$
</div>

# Convection Term: Dirichlet Boundary Condition
If a Dirichlet boundary condition is imposed on a face, the expression for $$\phi_k$$ is quite simple: it is the value that the boundary condition specifies. Be careful, though. Its value must be calculated at the midpoint of the face. Because $$\phi_k$$ is known, $$F^{\mathrm{c}}_k$$ does not contain any unknown and is readily available. Therefore, it can be moved to $$\textbf{b}$$.

<div class="post-blue">
The convection flux for Dirichlet faces is thus:
$$
F^{\mathrm{c}}_k = \rho \phi_k \textbf{v}_k \cdot \textbf{S}_k
$$

and since all the parameters are known, it does not have any contribution to \(\textbf{A}\). This known flux is moved to the right-hand side vector:

$$
\textbf{b}(i) \mathrel{+}= -\rho \phi_k \textbf{v}_k \cdot \textbf{S}_k
$$
</div>

# Convection Term: Neumann Boundary Condition
In our case, $$\frac{\partial \phi}{\partial x} = 0$$ (outlet) is imposed on the right boundary, and $$\frac{\partial \phi}{\partial y} = 0$$ (y-symmetry) is imposed on the bottom boundary. For either case, given that the mesh is fine enough, we can simply say that the face-centered value is approximately the cell-centered value of the boundary cell, i.e.,

$$
\phi_k \approx \phi_{\mathrm{P}}
$$

<div class="post-blue">
The convection flux for Neumann faces is thus:
$$
F^{\mathrm{c}}_k = \rho \phi_k \textbf{v}_k \cdot \textbf{S}_k \approx \rho \phi_{\mathrm{P}} \textbf{v}_k \cdot \textbf{S}_k
$$

and it only has one contribution to the linear system:

$$
\textbf{A}(i,i) \mathrel{+}= \rho\textbf{v}_k \cdot \textbf{S}_k
$$
</div>

## Diffusion Term
Denote the diffusion term of each face $$k$$ as:

$$ F^{\mathrm{d}}_k = \Gamma (\nabla \phi)_k \cdot \textbf{S}_k $$

The dot product between the gradient and the face normal is the directional derivative in the direction of the face normal:

$$
F^{\mathrm{d}}_k = \Gamma (\nabla \phi)_k \cdot \textbf{S}_k = \Gamma \{(\nabla \phi)_k \cdot \textbf{n}\} S_k = \Gamma \frac{\partial\phi_k}{\partial \textbf{n}} S_k
$$

$$\Gamma$$ and $$S_k$$ are known, so what's left is to approximate the directional derivative at the midpoint of face $$k$$ with cell centered values.

# Diffusion Term: Internal Faces
If the face midpoint lies on the segment that joins the straddling cell centers, we can use the central difference:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} \approx \frac{\phi_{\mathrm{N}_k}-\phi_{\mathrm{P}}}{|\bm{r}_{\mathrm{N}_k} - \bm{r}_{\mathrm{P}}|}
$$

However, as in the case of the convection term, the face midpoint might not lie on the segment. Two auxiliary points, $$\mathrm{P'}$$ and $$\mathrm{N}'_{k}$$ are introduced:

|![Diffusion term internal face auxiliary points](/images/transport-equation-fvm-unstructured-mesh/fvm_cell_diffustion_aux.svg)
|:--:| 
| *Figure 3:* The auxiliary points $$\mathrm{P'}$$ and $$\mathrm{N}'_{k}$$ |

These points lie on the line that is parallel to the face normal $$\textbf{n}$$ and intersects the face mid-point. Therefore, we have this central difference:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} \approx \frac{\phi_{\mathrm{N}'_k}-\phi_{\mathrm{P'}}}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

The scalar value at these auxiliary points is expressed as:

$$
\phi_{\mathrm{P'}} \approx \phi_{\mathrm{P}} + (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})
$$

$$
\phi_{\mathrm{N}'_k} \approx \phi_{\mathrm{N}_k} + (\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k})
$$

The central difference expressed with cell-centered values is therefore:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} \approx \frac{\phi_{\mathrm{N}_k}-\phi_{\mathrm{P}}}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|} + \frac{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

The second term is again treated as deferred correction.

The position of the auxiliary points are obtained by projecting $$\bm{r}_k - \bm{r}_{\mathrm{P}}$$ and $$\bm{r}_{\mathrm{N}_k} - \bm{r}_k$$ onto the face normal $$\textbf{n}$$. The magnitude of the shorter projection is recorded as $$\beta$$:

$$
\beta = \mathrm{min}\left((\bm{r}_k - \bm{r}_{\mathrm{P}}) \cdot \textbf{n}, (\bm{r}_{\mathrm{N}_k} - \bm{r}_k) \cdot \textbf{n}\right)
$$

and we use this magnitude to scale $$\textbf{n}$$ into cell $$\mathrm{P}$$ to get $$\textbf{r}_{\mathrm{P'}}$$ and into cell $$\mathrm{N}'_k$$ to get $$\textbf{r}_{\mathrm{N}'_k}$$:

$$
\textbf{r}_{\mathrm{P'}} = \textbf{r}_k - \beta\textbf{n}
$$

$$
\textbf{r}_{\mathrm{N}'_k} = \textbf{r}_k + \beta\textbf{n}
$$

<div class="post-blue">
The diffusion flux for internal faces is thus:

$$
F^{\mathrm{d}}_k = \Gamma \left\{\frac{\phi_{\mathrm{N}_k}-\phi_{\mathrm{P}}}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|} + \frac{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}\right\} S_k
$$

Collecting the coefficients:
$$
F^{\mathrm{d}}_k = -\left\{\frac{\Gamma S_k}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}\right\}\phi_{\mathrm{P}} + \left\{\frac{\Gamma S_k}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}\right\}\phi_{\mathrm{N}_k} + \frac{\Gamma S_k\{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

Its contribution to the linear system:
$$
\textbf{A}(i,i) \mathrel{+}= \frac{\Gamma S_k}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

$$
\textbf{A}(i,j) \mathrel{+}= -\frac{\Gamma S_k}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

$$
\textbf{b}(i) \mathrel{+}= \frac{\Gamma S_k\{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}}{|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|}
$$

<i>Note that the signs of the contribution from the diffusion terms are flipped because the total diffusion flux is being subtracted from the total convection flux in the transport equation.</i>
</div>

# Diffusion Term: Dirichlet Boundary Condition
Implementing a Dirichlet boundary condition for the diffusion term starts with the one-sided difference:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} \approx \frac{\phi_k-\phi_{\mathrm{P}}}{|\bm{r}_k - \bm{r}_{\mathrm{P}}|}
$$

Once again, $$\bm{r}_k - \bm{r}_{\mathrm{P}}$$ is likely not parallel to $$\textbf{n}$$:

|![Diffusion term Dirichlet boundary auxiliary points](/images/transport-equation-fvm-unstructured-mesh/fvm_cell_diffusion_dirichlet.svg)
|:--:| 
| *Figure 4:* Auxiliary point correction for $$\frac{\partial\phi_k}{\partial \textbf{n}}$$ |

The same auxiliary-point idea can be used to correct this:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} \approx \frac{\phi_k-\phi_{\mathrm{P'}}}{|\bm{r}_k - \bm{r}_{\mathrm{P}'}|} = \frac{\phi_k-\phi_{\mathrm{P}}}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|} - \frac{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}
$$

where $$\phi_k$$ is known from the boundary condition and the last term is a deferred correction. $$\textbf{r}_{\mathrm{P'}}$$ is obtained by projecting $$\textbf{r}_k - \textbf{r}_{\mathrm{P}}$$ onto $$\textbf{n}$$:

$$
\textbf{r}_{\mathrm{P'}} = \textbf{r}_k - \{(\textbf{r}_k - \textbf{r}_{\mathrm{P}}) \cdot \textbf{n}\}\textbf{n}
$$

<div class="post-blue">
The diffusion flux for Dirichlet faces is thus:
$$
F^{\mathrm{d}}_k = \Gamma \left\{\frac{\phi_k-\phi_{\mathrm{P}}}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|} - \frac{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}\right\} S_k
$$

Collecting the coefficients:
$$
F^{\mathrm{d}}_k = -\left\{\frac{\Gamma S_k}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}\right\}\phi_{\mathrm{P}} + \left\{\frac{\Gamma S_k}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}\right\}\phi_k - \frac{\Gamma S_k\{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}
$$

Its contribution to the linear system:
$$
\textbf{A}(i,i) \mathrel{+}= \frac{\Gamma S_k}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}
$$

$$
\textbf{b}(i) \mathrel{+}= \left\{\frac{\Gamma S_k}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}\right\}\phi_k - \frac{\Gamma S_k\{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}}{|\bm{r}_k - \bm{r}_{\mathrm{P'}}|}
$$
</div>

# Diffusion Term: Neumann Boundary Condition
Because the two Neumann boundary conditions mentioned earlier are applied along the right ($$\frac{\partial \phi}{\partial x} = 0$$) and the bottom ($$\frac{\partial \phi}{\partial y} = 0$$) boundary of the domain, the fact that the face normal $$\textbf{n}$$ is directed along the x-axis for the right boundary and along the y-axis for the bottom boundary allows us to say:

$$
\frac{\partial\phi_k}{\partial \textbf{n}} = \frac{\partial \phi_k}{\partial \hat{\textbf{x}}} = \frac{\partial \phi_k}{\partial x} = 0 \text{ (right)}
$$

and

$$
\frac{\partial\phi_k}{\partial \textbf{n}} = -\frac{\partial \phi_k}{\partial \hat{\textbf{y}}} = -\frac{\partial \phi_k}{\partial y} = 0 \text{ (bottom)}
$$

<div class="post-blue">
The diffusion flux for Neumann faces is thus:
$$
F^{\mathrm{d}}_k = 0
$$
so it does not contribute anything to the linear system.
</div>

[rans]: https://en.wikipedia.org/wiki/Reynolds-averaged_Navier%E2%80%93Stokes_equations
[ferziger]: https://link.springer.com/book/10.1007/978-3-642-56026-2
[moukalled]: https://link.springer.com/book/10.1007/978-3-319-16874-6