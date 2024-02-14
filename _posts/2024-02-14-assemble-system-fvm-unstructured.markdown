---
layout: post
title:  "Assembling Linear System from Unstructured Finite-Volume Discretization"
date:   2024-02-14 11:50:00 +0700
categories: cfd
---
## Why?
This note presents an algorithm for assembling the linear system that is a result of the discretization of the scalar transport equation on an unstructured mesh. I wrote this while researching the finite-volume method. The discretization process was discussed in the [previous post][fvm_post].

## Data Structures
First, we need mesh information. The structures that hold this information is really simple. I store the nodes, the faces, and the cells (we will collectively call these *entities*) in their own arrays (`std::vector` in my case). They are also ordered because I need to access each entity by its index in the array.

I am currently relying on [Gmsh][gmsh] to generate finite-volume mesh for me. Sadly, since Gmsh is a finite-element mesh generator, it does not support finite-volume mesh by default, so connectivity information which is important for finite-volume algorithms is not exported from the program. I wrote this [simple Python script][convert_gmsh_fvm] to extract this information from the original finite-element mesh file and store it in three files: `node`, `face`, and `cell`. I then just read these files and store the data in the three vectors.

Nodes, faces, and cells are represented as objects, each with their own information. A `Node` only carries its coordinates in the Cartesian grid. A `Face` has the nodes it connects and the cells it belongs to. In a 2D triangular mesh, only two cells straddle an internal face. One cell is regarded as the *owner* of the face and the other the *neighbor*. What makes a cell the owner cell is arbitrary and is not important. What matters is that there should be only one owner cell and that a face vector *always* points away from the owner cell. If this very important property is neglected, the sign of some convection fluxes would be incorrect!

Finally, a `Cell` has its constituent nodes and faces. For each face, `Cell` also has the index of the neighbor cell associated to that face. Both `Face` and `Cell` also carry their geometrical properties like centroid, length, area, and volume, depending on the dimensionality of the simulation.

## Linear System
I need this linear system:

$$
\left(a_{\mathrm{P}}\phi_{\mathrm{P}} + \sum_k a_{\mathrm{N}_k}\phi_{\mathrm{N}_k}\right)_i = \mathrm{RHS}_i
$$

which is the scalar transport equation written with a sum of fluxes in and out of each cell, $$F^{\mathrm{c}} - F^{\mathrm{d}} = 0$$. Constructing the linear system is a matter of correctly calculating the coefficients $$a$$ and the constant terms and putting them in the right places.

## Algorithm
Suppose that we are working with a 2D triangular mesh. In this mesh, an internal cell has exactly 3 neighbors. A boundary cell has at least one and at most two neighbor cells. We loop through the cell arrays and construct the equations cell-by-cell. The algorithm below assumes that $$\nabla \phi$$ is available. If it's not, we can first set it to zero. Also, the indices here are zero-based.

<div class="post-box green">
<b>Assembly algorithm:</b><br>
\(n \leftarrow\) number of cells<br>
Zero-initialize matrix \(\textbf{A}_{n \times n}\) and vector \(\textbf{b}_{n \times 1}\)<br>
For each cell \(i = 0 \rightarrow n-1\):<br>
<img src="/images/indent.svg">For each of its faces \(k = 0 \rightarrow 3\):<br>
<img src="/images/indent2.svg">If it's an internal face:<br>
<img src="/images/indent3.svg"><span class="comment">Determine the neighbor cell index</span><br>
<img src="/images/indent3.svg">\(j \leftarrow\) cell neighbor index associated to face \(k\)<br>
<img src="/images/indent3.svg"><span class="comment">Store convection contribution</span><br>
<img src="/images/indent3.svg">\(\textbf{S}_k \leftarrow\) the vector of face \(k\)<br>
<img src="/images/indent3.svg">If cell \(i\) is not the owner of face \(j\), flip \(\textbf{S}_k\)<br>
<img src="/images/indent3.svg">Store \(\textbf{v}_k \cdot \textbf{S}_k\) to avoid repetitive computation<br>
<img src="/images/indent3.svg">\(\textbf{A}(i,i) \mathrel{+}= \rho(1-\xi_k)(\textbf{v}_k \cdot \textbf{S}_k)\)<br>
<img src="/images/indent3.svg">\(\textbf{A}(i,j) \mathrel{+}= \rho\xi_k(\textbf{v}_k \cdot \textbf{S}_k)\)<br>
<img src="/images/indent3.svg">\(\textbf{b}(i) \mathrel{+}= -\rho(\textbf{v}_k \cdot \textbf{S}_k)(\nabla \phi)_{k'} \cdot (\bm{r}_k-\bm{r}_{k'})\)<br>
<img src="/images/indent3.svg"><span class="comment">Store diffusion contribution</span><br>
<img src="/images/indent3.svg">\(a^d_\mathrm{P} \leftarrow \Gamma S_k/|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|\)<br>
<img src="/images/indent3.svg">\(\textbf{A}(i,i) \mathrel{+}= a^d_\mathrm{P}\)<br>
<img src="/images/indent3.svg">\(\textbf{A}(i,j) \mathrel{+}= -a^d_\mathrm{P}\)<br>
<img src="/images/indent3.svg">\(\textbf{b}(i) \mathrel{+}= \Gamma S_k\{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}/|\bm{r}_{\mathrm{N}'_k} - \bm{r}_{\mathrm{P'}}|\)<br>
<img src="/images/indent3.svg"><span class="comment">Done, continue face loop</span><br>
<img src="/images/indent2.svg">If it's a boundary face:<br>
<img src="/images/indent3.svg">\(\textbf{S}_k \leftarrow\) the vector of face \(k\)<br>
<img src="/images/indent3.svg">No flipping required since there is only the owner cell<br>
<img src="/images/indent3.svg">If boundary condition type is Dirichlet:<br>
<img src="/images/indent4.svg"><span class="comment">Store convection contribution</span><br>
<img src="/images/indent4.svg">\(\textbf{b}(i) \mathrel{+}= -\rho \phi_k \textbf{v}_k \cdot \textbf{S}_k\)<br>
<img src="/images/indent4.svg"><span class="comment">Store diffusion contribution</span><br>
<img src="/images/indent4.svg">\(\textbf{A}(i,i) \mathrel{+}= \Gamma S_k/|\bm{r}_k - \bm{r}_{\mathrm{P'}}|\)<br>
<img src="/images/indent4.svg">\(\textbf{b}(i) \mathrel{+}= \left\{\Gamma S_k/|\bm{r}_k - \bm{r}_{\mathrm{P'}}|\right\}\phi_k\) \(- \Gamma S_k\{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}/|\bm{r}_k - \bm{r}_{\mathrm{P'}}|\)<br>
<img src="/images/indent4.svg"><span class="comment">Done, continue face loop</span><br>
<img src="/images/indent3.svg">If boundary condition type is Neumann:<br>
<img src="/images/indent4.svg"><span class="comment">Store convection contribution</span><br>
<img src="/images/indent4.svg">\(\textbf{A}(i,i) \mathrel{+}= \rho\textbf{v}_k \cdot \textbf{S}_k\)<br>
<img src="/images/indent4.svg"><span class="comment">No diffusion contribution</span><br>
<img src="/images/indent4.svg"><span class="comment">Done, continue face loop</span><br>
</div>

[fvm_post]: {% post_url 2024-02-08-transport-equation-fvm-unstructured-mesh %}
[gmsh]: https://gmsh.info/
[convert_gmsh_fvm]: https://github.com/lum-supakorn/convert_gmsh_fvm
