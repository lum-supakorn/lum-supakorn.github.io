---
layout: post
title:  "Gauss-Seidel as a Projection Method"
date:   2024-02-22 11:20:00 +0700
categories: numerics
---

I learned from Yousef Saad's [Iterative Methods for Sparse Linear Systems](https://www-users.cse.umn.edu/~saad/books.html) that the Gauss-Seidel method is just a projection method where the subspace is $$\mathrm{span}\{\textbf{e}_i\}$$. This was surprising to me because the Gauss-Seidel algorithm I learned from my introductory numerical methods class did not seem to involve any projection. I wanted to see how this is the case and this note is the result of this little exploration.

A general projection method for solving a linear system iteratively is defined as:

<span class="center">
Find $$\tilde{\textbf{x}} \in \textbf{x}_0 + \mathcal{K}$$,<img src="/images/indent.svg"/>such that $$\textbf{b}-\textbf{A}\tilde{\textbf{x}} \perp \mathcal{L}$$
</span>

where $$\tilde{\textbf{x}}$$ is the approximate solution we are looking for, $$\textbf{x}_0$$ is an initial guess, $$\mathcal{K}$$ is the search subspace, $$\textbf{A}$$ and $$\textbf{b}$$ define our linear system, and $$\mathcal{L}$$ is the left subspace. The distinction between $$\mathcal{K}$$ and $$\mathcal{L}$$ does not matter in Gauss-Seidel case because they are both $$\mathrm{span}\{\textbf{e}_i\}$$. In any iteration $$k$$ we are looking for $$\textbf{x}_{k+1}$$ in the affine space $$\textbf{x}_k + \mathrm{span}\{\textbf{e}_i\}$$ and $$i$$ is cycled through $$1,\dots,n$$ where $$n$$ is the dimension of the problem.

Let's start with $$i = 1$$. In the first iteration we are looking for $$\textbf{x}_1$$ in $$\textbf{x}_0 + \mathrm{span}\{\textbf{e}_1\}$$ (line through $$\textbf{x}_0$$ along the direction of $$\textbf{e}_1$$) such that $$\textbf{r}_1 = \textbf{b}-\textbf{A}\textbf{x}_1 \perp \mathrm{span}\{\textbf{e}_1\}$$. That gives us:

$$
\begin{aligned}
(\textbf{b}-\textbf{A}\textbf{x}_1)^\top \textbf{e}_1 &= 0 \\
\textbf{b}^\top \textbf{e}_1-(\textbf{A}\textbf{x}_1)^\top \textbf{e}_1 &= 0 \\
b_1-\textbf{x}_1^\top \textbf{A}^\top \textbf{e}_1 &= 0 \\
\textbf{x}_1^\top \textbf{A}^\top \textbf{e}_1 &= b_1 \\
\textbf{x}_1^\top
\begin{bmatrix}
a_{11} \\
a_{12} \\
\vdots \\
a_{1n}
\end{bmatrix} &= b_1 \\
a_{11} x_{1,1} + a_{12} x_{1,2} + \dots + a_{1n} x_{1,n} &= b_1
\end{aligned}
$$

Since we are looking for $$\textbf{x}_1$$ strictly along $$\textbf{x}_0 + \mathrm{span}\{\textbf{e}_1\}$$, $$\textbf{x}_1$$ differs from $$\textbf{x}_0$$ only in the $$\textbf{e}_1$$-component. In other words, only $$x_{1,1}$$ must be looked for, and rest of the components are from the previous iteration. We now have:

$$
x_{1,1} = \frac{b_1 - \sum_{j=2}^{n}a_{1j} x_{0,j}}{a_{11}}
$$

After $$\textbf{x}_1$$ is found, $$\textbf{x}_2$$ can be searched for along $$\textbf{x}_1 + \mathrm{span}\{\textbf{e}_2\}$$. This process goes on until convergence is reached. In iteration $$k$$ we update the $$i$$th component of $$\tilde{\textbf{x}}$$ with this expression:

$$
x_{k,i} = \frac{b_i - \sum_{\substack{j=1\\j\neq i}}^{n}a_{ij} x_{k-1,j}}{a_{ii}}
$$

Now this starts to look like the Gauss-Seidel algorithm that I was taught. I implemented this on a 2-by-2 linear system to visualize the search:

![Gauss-Seidel Projection Search](/images/gauss-seidel-projection/gs-1.png)

And of course more iterations yield a more accurate result:

![Gauss-Seidel Projection Search](/images/gauss-seidel-projection/gs-2.png)