This folder contains the first steps in the development of magneto-fluid code. Here we are only working with the fluid component in 2D. We include the `*.cpp` files that should be complied to get your own `*.mex` files that can be used by MATLAB. The two approaches to modelling 2D Euler equations are:
- `LaxFriedrich_2DFlow2.cpp` routine solves the 2D Euler equations using the Lax-Friedrich scheme.
- `MacCormack_2DFlow8.cpp` routine solves the 2D Euler equations using the MacCormack scheme.

The docstrings of the scripts provide the documentation which explains how to run them.