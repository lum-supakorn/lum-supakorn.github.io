---
layout: post
title:  "Strassen's Matrix Multiplication Algorithm"
date:   2024-02-25 15:25:00 +0700
categories: numerics
---

This note is based on problem 32.2 of Lloyd N. Trefethen and David Bau III's [Numerical Linear Algebra](https://www.amazon.com/Numerical-Linear-Algebra-Lloyd-Trefethen/dp/0898713617).

Strassen's algorithm is a matrix multiplication algorithm that is faster than the usual approach. Suppose we have this $$m \times m$$ square matrix multiplication:

$$
\begin{bmatrix}
W & X \\
Y & Z
\end{bmatrix} =
\begin{bmatrix}
A & B \\
C & D
\end{bmatrix}
\begin{bmatrix}
E & F \\
G & H
\end{bmatrix}
$$

where the matrices are partitioned into square blocks of size $$\frac{m}{2} \times \frac{m}{2}$$. This partitioning is [conforming](https://en.wikipedia.org/wiki/Block_matrix), so we can carry out the multiplication as if the blocks are matrix elements:

$$
\begin{bmatrix}
W & X \\
Y & Z
\end{bmatrix} =
\begin{bmatrix}
AE+BG & AF+BH \\
CE+DG & CF+DH
\end{bmatrix}
$$

This usual approach requires 4 matrix additions and 8 matrix multiplications. For an $$m \times m$$ matrix, addition/subtraction requires $$m^2$$ operations, and multiplication requires $$2m^3-m^2$$ operations. So in total the usual multiplication requires $$16m^3 - 4m^2$$ or $$O(m^3)$$ operations.

[Volker Strassen](https://en.wikipedia.org/wiki/Volker_Strassen) came up with the following multiplication process:

$$
\begin{aligned}
    P_1 &= (A+D)(E+H)\\
    P_2 &= (C+D)E\\
    P_3 &= A(F-H)\\
    P_4 &= D(G-E)\\
    P_5 &= (A+B)H\\
    P_6 &= (C-A)(E+F)\\
    P_7 &= (B-D)(G+H)\\
    W &= P_1 + P_4 - P_5 + P_7\\
    X &= P_3 + P_5\\
    Y &= P_2 + P_4\\
    Z &= P_1 + P_3 - P_2 + P_6
\end{aligned}
$$

This process requires 18 matrix additions/subtractions but only 7 multiplications. If the multiplications to get all the $$P$$ matrices are carried out conventionally, this process is still $$O(m^3)$$. However, if we also use Strassen's multiplication for the $$P$$'s we end up with a recursive process that shows better performance.

The process of multiplying two $$m \times m$$ matrices, $$M_1$$ and $$M_2$$, where $$m = 2^k, k = 1,2,3,\dots$$ is as follows:
<div class="post-box green">
<b>Recursive Strassen's Matrix Multiplication Algorithm:</b><br>
function \(\text{strassen}(M_1,M_2)\):<br>
<span class="indent-1">\(m \leftarrow\) size of matrix  \(M_1\) and \(M_2\)<br></span>
<span class="indent-1">if \(m = 2\):<br></span>
<div class="algo-box indent-2">
    <span class="comment">Base case, perform normal matrix multiplication</span><br>
    return \(M_1M_2\)
</div>
<span class="indent-1">else:<br></span>
<div class="algo-box indent-2">
    <span class="comment">Recursive case, perform Strassen's matrix multiplication</span><br>
    Partition \(M_1\) and \(M_2\) into \(\frac{m}{2}\times\frac{m}{2}\) blocks:<br>
    \(M_1 = \begin{bmatrix}
    A & B \\
    C & D
    \end{bmatrix}\) and \(M_2 = \begin{bmatrix}
    E & F \\
    G & H
    \end{bmatrix}\)<br>
    <span class="comment">Follow Strassen's process</span><br>
    \(
    \begin{aligned}
        P_1 &\leftarrow \text{strassen}(A+D,E+H)\\
        P_2 &\leftarrow \text{strassen}(C+D,E)\\
        P_3 &\leftarrow \text{strassen}(A,F-H)\\
        P_4 &\leftarrow \text{strassen}(D,G-E)\\
        P_5 &\leftarrow \text{strassen}(A+B,H)\\
        P_6 &\leftarrow \text{strassen}(C-A,E+F)\\
        P_7 &\leftarrow \text{strassen}(B-D,G+H)\\
        W &\leftarrow P_1 + P_4 - P_5 + P_7\\
        X &\leftarrow P_3 + P_5\\
        Y &\leftarrow P_2 + P_4\\
        Z &\leftarrow P_1 + P_3 - P_2 + P_6
    \end{aligned}
    \)<br>
    return \(\begin{bmatrix}
    W & X \\
    Y & Z
    \end{bmatrix}
    \)
</div>
</div>

Following this recursive process shows us that it requires

$$12\cdot7^{(k-1)} + 18\sum_{i=1}^{k-1}7^{i-1}2^{2(k-i)}$$

operations. The leading term is $$7^{(k-1)}$$ so this algorithm is $$O(7^k)$$ or, through [logarithm power change](https://math.stackexchange.com/questions/453918/proof-of-logarithm-power-change), $$O(m^{\log_2{7}}) \approx O(m^{2.81})$$ which is better than $$O(m^3)$$.