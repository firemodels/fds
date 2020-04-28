#!/bin/bash

# The -O option runs cases
#  casea.fds, ..., caseh.fds all on the same node using
#  1, ..., 8 threads respectively where case is
#  openmp_test64 and openmp_test128

$QFDS -t -O 8 -A -d Timing_Benchmarks openmp_test64

$QFDS -t -O 8 -A -d Timing_Benchmarks openmp_test128
