#!/bin/bash

# The -O option runs cases
#  casea.fds, ..., caseh.fds all on the same node using
#  1, ..., 8 threads respectively where case is
#  openmp_test64 and openmp_test128

$QFDS -t -o 1 -A -d Timing_Benchmarks openmp_test64a.fds
$QFDS -t -o 2 -A -d Timing_Benchmarks openmp_test64b.fds
$QFDS -t -o 3 -A -d Timing_Benchmarks openmp_test64c.fds
$QFDS -t -o 4 -A -d Timing_Benchmarks openmp_test64d.fds
$QFDS -t -o 5 -A -d Timing_Benchmarks openmp_test64e.fds
$QFDS -t -o 6 -A -d Timing_Benchmarks openmp_test64f.fds
$QFDS -t -o 7 -A -d Timing_Benchmarks openmp_test64g.fds
$QFDS -t -o 8 -A -d Timing_Benchmarks openmp_test64h.fds

$QFDS -t -o 1 -A -d Timing_Benchmarks openmp_test128a.fds
$QFDS -t -o 2 -A -d Timing_Benchmarks openmp_test128b.fds
$QFDS -t -o 3 -A -d Timing_Benchmarks openmp_test128c.fds
$QFDS -t -o 4 -A -d Timing_Benchmarks openmp_test128d.fds
$QFDS -t -o 5 -A -d Timing_Benchmarks openmp_test128e.fds
$QFDS -t -o 6 -A -d Timing_Benchmarks openmp_test128f.fds
$QFDS -t -o 7 -A -d Timing_Benchmarks openmp_test128g.fds
$QFDS -t -o 8 -A -d Timing_Benchmarks openmp_test128h.fds
