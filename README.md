# General

A re-write of an old, badly written, program using templates.  I plan to use this program
to practice C++20 concepts, constraints, and ranges.

This quick program implements a Runge-Kutta-Fhelberg 
method with a rel_err=0.  Essentially this is the GSL method gsl_odeiv2_step_rkf45
with rel_err = 0.

Note, I still need to greatly refactor and optimize this program because it's mainly written to practice C++20 and use concepts (which is also in progress).

# Usage

The program solves the system:

dy_i(x_j)/dt = F_i(t, x_j, y_j)

see the example for how this is implemented.

## Compilation

To run the test enter the command:

g++ -O3 -std=c++11 newton_test.cc

To verify the output should be:

-1.49956e+11,2.10468e+10,0,-4531.51,-29292.7,0,


