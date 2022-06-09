## Thread Checker

This folder contains two simple, identical FDS input files named `race_test_1.fds` and `race_test_4.fds`. The former is run with one core using the non-OpenMP version of FDS. The latter is run with four cores using the OpenMP version. The metric of success is an identical time history of both temperature and velocity at a given point in the domain.

If this simple test fails, you will need to use Intel Inspector to find out where in the code there is a race condition. A brief tutorial on Inspector can be found [here](https://github.com/firemodels/fds/wiki/Intel-Diagnostic-Tools).

