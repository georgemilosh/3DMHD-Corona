# 3DMHD-Corona
3D MHD model of solar coronal loop structure that was developed for a ![master thesis](/docs/thesis.pdf) by [George Miloshevich](https://georgemilosh.github.io) that was defended in summer 2012.

![Upflow trapped in magnetic arcade](/docs/out.gif)

## Introduction
This thesis explores thermalization of chromospheric upflows as they reach solar corona and are trapped in the magnetic field loops. The crucial mechanism behind this transfer of kinetic energy to heat is the ability of supersonic flows to steepen and generate shorter scales which can dissipate. The full treatment would require at least two-fluid effects which have the relevant small scales. However for simulation purposes we will limit ourselves to 3D MHD model of supersonic flows with adiabatic isothermal closure and no Hall term and leave these effects for future developments. Then the primary challenge is to accurately model such high Mach number flows which we achieve by using a combination of implicit MacCormack scheme and artificial viscosity. The project culminates in the illustration of coronal loop structure formation and heating via numerical simulation. 


## Technical Details

The simulation is written in C embedded in MATLAB 2012 for flexibility. C library uses `mex.h` which is provided via MATLAB. For more information about how to use this approach see [mex files](https://www.mathworks.com/help/matlab/call-mex-file-functions.html). The point to retain is that the simulation can be launched directly from MATLAB as though it was supported by MATLAB and the output is written to disk and RAM. 

## Numerics

The scripts rely on MacCormack method for time integration and artificial viscosity for shock capturing. The adaptation is based on the famous book [Numerical Recipes in C](https://dl.acm.org/doi/10.5555/148286). The MacCormack method is a two-step predictor-corrector method that is second order accurate in both space and time. The artificial viscosity is a numerical technique that is used to capture shocks. The code runs on a single core. 