# General

A re-write of an old, program using templates and concepts.  I plan to use this program to practice C++20 concepts, coroutines, and ranges. The plan is to create a general solver class which takes as input a Butcher Tableau dependency (as opposed to the current fixed constant routine).

The current program implements a Runge-Kutta-Fhelberg method with a rel_err=0.  Essentially this is the GSL method gsl_odeiv2_step_rkf45
with rel_err = 0.


# Usage

See `newton_test.cpp` for a basic example. The unit tests also have variations on the example. The RK45 class solves the system:

```
dy_i(x_j)/dt = F_i(t, x_j, y_j)
```

using the [Runge–Kutta–Fehlberg method](https://en.wikipedia.org/wiki/Runge%E2%80%93Kutta%E2%80%93Fehlberg_method).

# Compilation

## Cmake

Standard cmake commands

```
mkdir build
cd build
cmake ..
make
```

Whcih will also run tests. To examine the tests with more detail, run the binaries in the `tests` directory.


## CLI

The code only compiles for compilers with c++20 support.  If you have gcc-10, or greater, the command to compile is:

```
g++-10 -O3 -std=c++20 newton_test.cpp -I rk45
```

To verify the output should be run './a.out' and the last line should be:

```
Newton solution
-1.49956e+11,2.10468e+10,0,-4531.51,-29292.7,0,
Total runs: 5000
Average duration run: 70 microseconds
```

The current test uses a gcc output to test but I should use an analytically solvable model to test against at some point.



# Performance

Currently it's a factor of 1.05 times slower than GSL runtimes (see gcc example) on a intel i7.

