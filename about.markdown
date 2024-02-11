---
layout: page
title: About
permalink: /about/
---

![Lum](/images/lum.jpg)

Hi, my name is Lum. I'm an engineer from Thailand.\
I love programming to solve problems in maths and physics.\
I have to self-teach many of the things I want to learn and I find writing about them helps.

I'm currently working as a computational engineer at [HiveGround](https://www.hiveground.com/).

<br><br>

# Things I've worked on
- 2023
  - Designed and conducted an experiment to collect aerodynamics forces and moments on [VETAL][VETAL] during hover with a 6-DoF load sensor. Learned a lot about how to report measurement uncertainty mostly from [GUM 1995][GUM] and [John R. Taylor's error analysis book][error].
  - Made a Flask + Bootstrap web application for automatic aircraft component flight time tracking
  - Collected spatial slipstream velocity data induced by a drone propeller for performance evaluation and CFD validation.
  - Trained another intern! (more mechanical and electrical stuff for aerodynamics experiments).
  - Learned the basic mathematics beind finite-element methods from [MIT 18.085][strang] and [MIT 16.920J][pde]. Intended to apply the technique to a 3D incompressible code but abandoned the approach once I realized I needed accurate trailing edge separation physics.
  - Started digging into finite-volume incompressible flow algorithms.
- 2022
  - Wrote a Python pipeline that turns rough 3D-scan data of drone propellers into smooth (with convex hull algorithm) airfoil sections that were used to re-construct the propellers in OpenVSP into a smooth geometry model which can be meshed nicely with Pointwise.
  - Attempted to make an agricultural drone sprayer simulation by using particle dynamics with drone slipstream CFD data and numerical integration. Terminated because the simulation needed time-accurate CFD data which was prohibitive at the time.
  - Trained my first intern (mechanical and electrical stuff for aerodynamics experiments).
  - Completed the first iteration of propeller-wing incompressible flow simulation of VETAL with OpenFOAM and Pointwise.

[VETAL]: https://www.hiveground.com/vetal/
[GUM]: https://www.bipm.org/documents/20126/2071204/JCGM_100_2008_E.pdf/cb0ef43f-baa5-11cf-3f85-4dcd86f77bd6
[error]: https://www.amazon.com/Introduction-Error-Analysis-Uncertainties-Measurements/dp/093570275X
[strang]: https://ocw.mit.edu/courses/16-920j-numerical-methods-for-partial-differential-equations-sma-5212-spring-2003/
[pde]: https://ocw.mit.edu/courses/16-920j-numerical-methods-for-partial-differential-equations-sma-5212-spring-2003/