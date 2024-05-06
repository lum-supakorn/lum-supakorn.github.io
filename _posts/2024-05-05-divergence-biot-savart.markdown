---
layout: post
title:  "Divergence of Magnetic Field from the Biot-Savart Law and the Derivative of Definite Integrals"
date:   2024-05-05 23:00:00 +0700
categories: physics
---

### Divergence of Magnetic Field

I wanted to write about the fact that the divergence of the magnetic field is always zero (Gauss's law for magnetism):

$$
\nabla \cdot \textbf{B} = 0
$$

This means that the magnetic field is always divergence-free or *solenoidal*. Physically, this says that there is no *magnetic monopole*, because otherwise it would be a source or a sink of the magnetic field and the divergence would not be zero. Fair enough, but can we prove this with mathematics?

Griffiths' electrodynamics book explained this mathematically. It starts with the Biot-Savart law which describes the magnetic field induced by some medium with a steady volume current density:

$$
\textbf{B}(\textbf{r}) = \frac{\mu_0}{4\pi}\int_V{\frac{\textbf{J}\times\textbf{r}'}{|\textbf{r}'|^3}\mathrm{d}V}
$$

where $$\textbf{J}$$ is the current density of the infinitesimal volume and $$\textbf{r}'$$ is the displacement vector pointing from the the infinitesimal volume to the field point $$\textbf{r}$$. Now let's apply the divergence operator on this expression:

$$
\nabla \cdot \textbf{B}(\textbf{r}) = \nabla \cdot \left\{\frac{\mu_0}{4\pi}\int_V{\frac{\textbf{J}\times\textbf{r}'}{|\textbf{r}'|^3}\mathrm{d}V} \right\}
$$

Since the boundary of the volume integration is stationary, we can move the divergence operator inside the integral:

$$
\nabla \cdot \textbf{B}(\textbf{r}) = \frac{\mu_0}{4\pi}\int_V{\nabla \cdot \frac{\textbf{J}\times\textbf{r}'}{|\textbf{r}'|^3}\mathrm{d}V}
$$

We will now focus on the integrand. The displacement vectors are grouped together and the [divergence product rule](https://en.wikipedia.org/wiki/Divergence#Properties) is invoked:

$$
\begin{equation}
\nabla \cdot \left\{ \textbf{J}\times\frac{\textbf{r}'}{|\textbf{r}'|^3} \right\} = (\nabla \times \textbf{J}) \cdot \frac{\textbf{r}'}{|\textbf{r}'|^3} - \textbf{J} \cdot (\nabla \times \frac{\textbf{r}'}{|\textbf{r}'|^3})
\end{equation}
$$

Griffiths says that the first term is zero because $$\textbf{J}$$ only depends on the *primed* variables and invites us to solve Problem 1.63 to see that the second term is also zero (I will investigate this later). I am not satisfied with this explanation. What are *primed* and *unprimed* variables in the context of calculus? Does $$\textbf{J}$$ not vary in space throughout the medium so the spatial derivative of $$\textbf{J}$$ is non-zero? It seems like [a few people were confused too](https://www.reddit.com/r/AskPhysics/comments/gmke92/difference_between_primed_and_unprimed_variables/). This is what I would like to investigate.

### Derivative of Definite Integral

Let's say I have a function of $$x$$, $$f(x)$$. I ask: what is the derivative with respect to $$x$$ of a definite integral of $$f$$? If the boundary of integration $$[a,b]$$ gives us a constant definite integral value, the answer to that question is of course zero. That makes sense from this perspective, but let's see it from the Leibniz integral rule:

$$
\begin{equation}
\frac{\mathrm{d}}{\mathrm{d}x}\int_a^b{f(t)\mathrm{d}t} = \int_a^b{\frac{\mathrm{d}}{\mathrm{d}x}f(t)\mathrm{d}t} = \int_a^b{0\mathrm{d}t} = 0
\end{equation}
$$

This makes sense because I use $$t$$ as the variable of integration and the derivative with respect to $$x$$ yields zero. But what if I didn't do that? If I still use $$x$$ as the variable of integration:

$$
\begin{equation}
\frac{\mathrm{d}}{\mathrm{d}x}\int_a^b{f(x)\mathrm{d}x} = \int_a^b{\frac{\mathrm{d}}{\mathrm{d}x}f(x)\mathrm{d}x} = f(b)-f(a)
\end{equation}
$$

which is not necessarily zero! There must be a hole in my understanding of calculus here. The choice of this *variable of integration* should not affect the result. This confusion led me to [the definition of *Riemann integral*](https://en.wikipedia.org/wiki/Integral#Riemann_integral) which comes from the *Riemann sum*:

$$
\sum_{i=1}^n{f(t_i)\Delta_i}
$$

where the $$t_i$$'s are the *tags* of the sub-interval inside $$[a,b]$$. It became clear to me now that the value of the function $$f$$ inside the sum is the value of the function *at* the tag. It is locked in, so it becomes a constant and does not vary anymore. This is why it becomes zero in Equation 2. The writing of the definite integral in Equation 3 is probably vague or completely incorrect if we want to convey the derivative of the definite integral with constant boundary.

So the $$\textbf{J}$$ in Equation 1 is not a function that varies spatially in 3D space anymore, but a constant value that represents the volume current density of the *tagged* volume sub-interval in the Riemann sum. That is why the spatial derivative and therefore the curl of it is zero.