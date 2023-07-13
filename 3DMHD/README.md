Here the 3D MHD code is presented. We include `*.c` routines that should be complied to get your own `*.mex` files that can be used by MATLAB. 
- `MHD_MacCormack11.cpp` routine solves the 2D MHD equations using the MacCormack scheme.
- `MC3DMHD3.cpp` routines solves the 3D MHD equations using the MacCormack scheme.

The docstrings of the scripts provide the documentation which explains how to run them.
In addition we provide some MATLAB `*.m` files that are necessary to analyze the output of the simulation and prepare initial conditions.