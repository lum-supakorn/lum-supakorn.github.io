---
layout: post
title:  "Assembling Linear System from Unstructured Finite-Volume Discretization"
date:   2024-02-14 11:50:00 +0700
categories: cfd
---
## Why?
This note presents an algorithm for assembling the linear system that is a result of the discretization of the scalar transport equation on an unstructured mesh. I wrote this while researching the finite-volume method. The discretization process was discussed in the [previous post][fvm_post].

## Data Structures
First, we need mesh information. The structures that hold this information is really simple. I store the nodes, the faces, and the cells (we will collectively call these *entities*) in their own arrays. They are also ordered because I need to access each entity by its index in the array.

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

<input type="checkbox" id="geom" onclick="show('geom')" checked>
<label for="geom">Show geometry calculation</label>
<input type="checkbox" id="gradient" onclick="show('gradient')" checked>
<label for="gradient">Gradient is available</label>

<div class="post-box green">
<b>Assembly algorithm:</b><br>
\(n \leftarrow\) number of cells<br>
Zero-initialize matrix \(\textbf{A}_{n \times n}\) and vector \(\textbf{b}_{n \times 1}\)<br>
For each cell \(i = 0 \rightarrow n-1\):<br>
<img src="/images/indent.svg">For each of its faces \(k\):<br>
<img src="/images/indent2.svg">If it's an internal face:<br>
<div class="algo-box-group">
    <div class="algo-box indent-3">
    <span class="comment">Determine the neighbor cell index</span><br>
    \(j \leftarrow\) cell neighbor index associated to face \(k\)<br>
    </div>
    <div class="algo-box indent-3 geom">
    <span class="comment">Compute geometrical variables for convection</span><br>
    \(\textbf{S}_k \leftarrow\) the vector of face \(k\)<br>
    If cell \(i\) is not the owner of face \(k\), flip \(\textbf{S}_k\)<br>
    \(\textbf{v}_k \leftarrow\) velocity at face midpoint<br>
    Store \(\textbf{v}_k \cdot \textbf{S}_k\) to avoid repetitive computation<br>
    \(\textbf{r}_{k} \leftarrow\) the centroid of face \(k\)<br>
    \(\textbf{r}_\mathrm{P} \leftarrow\) the centroid of cell \(i\)<br>
    \(\textbf{r}_{\mathrm{N}_k} \leftarrow\) the centroid of cell \(j\)<br>
    \(\xi_k \leftarrow (\textbf{r}_k - \textbf{r}_{\mathrm{P}}) \cdot (\textbf{r}_{\mathrm{N}_k} - \textbf{r}_{\mathrm{P}})/|\textbf{r}_{\mathrm{N}_k} - \textbf{r}_{\mathrm{P}}|^2\)
    <span class="gradient"><br>\(\textbf{r}_{k'} \leftarrow (1-\xi_k)\textbf{r}_{\mathrm{P}} + \xi_k\textbf{r}_{\mathrm{N}_k}\)</span>
    </div>
    <div class="algo-box indent-3 gradient">
    <span class="comment">Interpolate gradients</span><br>
    \((\nabla \phi)_{k'} \leftarrow (1-\xi_k)(\nabla \phi)_\mathrm{P} + \xi_k(\nabla \phi)_{\mathrm{N}_k} \)<br>
    </div>
    <div class="algo-box indent-3">
    <span class="comment">Store convection contribution</span><br>
    \(\textbf{A}(i,i) \mathrel{+}= \rho(1-\xi_k)(\textbf{v}_k \cdot \textbf{S}_k)\)<br>
    \(\textbf{A}(i,j) \mathrel{+}= \rho\xi_k(\textbf{v}_k \cdot \textbf{S}_k)\)
    <span class="gradient"><br>\(\textbf{b}(i) \mathrel{+}= -\rho(\textbf{v}_k \cdot \textbf{S}_k)(\nabla \phi)_{k'} \cdot (\textbf{r}_k-\textbf{r}_{k'})\)</span><br>
    </div>
    <div class="algo-box indent-3 geom">
    <span class="comment">Compute geometrical variables for diffusion</span><br>
    \(\textbf{n} \leftarrow\) normalized \(\textbf{S}_k\)<br>
    \(\beta \leftarrow \mathrm{min}\left((\textbf{r}_k - \textbf{r}_{\mathrm{P}}) \cdot \textbf{n}, (\textbf{r}_{\mathrm{N}_k} - \textbf{r}_k) \cdot \textbf{n}\right)\)<br>
    \(\textbf{r}_{\mathrm{P'}} \leftarrow \textbf{r}_k - \beta\textbf{n}\)<br>
    \(\textbf{r}_{\mathrm{N}'_k} \leftarrow \textbf{r}_k + \beta\textbf{n}\)<br>
    </div>
    <div class="algo-box indent-3">
    <span class="comment">Store diffusion contribution</span><br>
    \(a^d_\mathrm{P} \leftarrow \Gamma S_k/|\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{P'}}|\)<br>
    \(\textbf{A}(i,i) \mathrel{+}= a^d_\mathrm{P}\)<br>
    \(\textbf{A}(i,j) \mathrel{+}= -a^d_\mathrm{P}\)
    <span class="gradient"><br>\(\textbf{b}(i) \mathrel{+}= \Gamma S_k\{(\nabla \phi)_{\mathrm{N}_k} \cdot (\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{N}_k}) - (\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}/|\textbf{r}_{\mathrm{N}'_k} - \textbf{r}_{\mathrm{P'}}|\)<br></span>
    </div>
</div>
<img src="/images/indent3.svg"><span class="comment">Done, continue face loop</span><br>
<img src="/images/indent2.svg">If it's a boundary face:<br>
<div class="algo-box-group">
    <div class="algo-box indent-3 geom">
    \(\textbf{S}_k \leftarrow\) the vector of face \(k\)<br>
    No flipping required since there is only the owner cell<br>
    \(\textbf{v}_k \leftarrow\) velocity at face midpoint<br>
    </div>
</div>
<img src="/images/indent3.svg">If boundary condition type is Dirichlet:<br>
<div class="algo-box-group">
    <div class="algo-box indent-4 geom">
    \(\textbf{r}_{k} \leftarrow\) the centroid of face \(k\)<br>
    \(\textbf{r}_\mathrm{P} \leftarrow\) the centroid of cell \(i\)<br>
    </div>
    <div class="algo-box indent-4">
    <span class="comment">Store convection contribution</span><br>
    \(\phi_k \leftarrow\) scalar value from boundary condition<br>
    \(\textbf{b}(i) \mathrel{+}= -\rho \phi_k \textbf{v}_k \cdot \textbf{S}_k\)<br>
    </div>
    <div class="algo-box indent-4">
    <span class="comment">Store diffusion contribution</span><br>
    \(\textbf{A}(i,i) \mathrel{+}= \Gamma S_k/|\textbf{r}_k - \textbf{r}_{\mathrm{P'}}|\)<br>
    \(\textbf{b}(i) \mathrel{+}= \left\{\Gamma S_k/|\textbf{r}_k - \textbf{r}_{\mathrm{P'}}|\right\}\phi_k\)
    <span class="gradient"><br><img src="/images/indent3.svg">\(- \Gamma S_k\{(\nabla \phi)_{\mathrm{P}} \cdot (\textbf{r}_{\mathrm{P'}} - \textbf{r}_{\mathrm{P}})\}/|\textbf{r}_k - \textbf{r}_{\mathrm{P'}}|\)</span>
    </div>
</div>
<img src="/images/indent4.svg"><span class="comment">Done, continue face loop</span><br>
<img src="/images/indent3.svg">If boundary condition type is Neumann:<br>
<div class="algo-box-group">
    <div class="algo-box indent-4">
    <span class="comment">Store convection contribution</span><br>
    \(\textbf{A}(i,i) \mathrel{+}= \rho\textbf{v}_k \cdot \textbf{S}_k\)<br>
    </div>
    <div class="algo-box indent-4">
    <span class="comment">No diffusion contribution</span><br>
    </div>
</div>
</div>

<script>
    function show(group) {
        let box = document.getElementById(group);
        let vis;
        if (box.checked) {
            vis = "inline";
        } else {
            vis = "none"
        }
        items = document.getElementsByClassName(group);
        for (let i = 0; i < items.length; i++) {
            items[i].style.display = vis;
        }
    }
</script>

[fvm_post]: {% post_url 2024-02-08-transport-equation-fvm-unstructured-mesh %}
[gmsh]: https://gmsh.info/
[convert_gmsh_fvm]: https://github.com/lum-supakorn/convert_gmsh_fvm
