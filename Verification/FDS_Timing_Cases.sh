#!/bin/bash
# cpu times are summed for all cases below and report in ... .
# cases should be representative of the kind of work
# performed by FDS.  It also should be likely that these cases 
# do not change. The idea is to reduce the peformance of FDS
# down to one number - a Dow Jones average so to speak
# These case below are a straw man
 
$QFDS -t -o 1 -d Timing_Benchmarks openmp_test64a.fds
$QFDS -t -o 2 -d Timing_Benchmarks openmp_test64b.fds
$QFDS -t -o 3 -d Timing_Benchmarks openmp_test64c.fds
$QFDS -t -o 4 -d Timing_Benchmarks openmp_test64d.fds
$QFDS -t -o 5 -d Timing_Benchmarks openmp_test64e.fds
$QFDS -t -o 6 -d Timing_Benchmarks openmp_test64f.fds
$QFDS -t -o 7 -d Timing_Benchmarks openmp_test64g.fds
$QFDS -t -o 8 -d Timing_Benchmarks openmp_test64h.fds

$QFDS -t -o 1 -d Timing_Benchmarks openmp_test128a.fds
$QFDS -t -o 2 -d Timing_Benchmarks openmp_test128b.fds
$QFDS -t -o 3 -d Timing_Benchmarks openmp_test128c.fds
$QFDS -t -o 4 -d Timing_Benchmarks openmp_test128d.fds
$QFDS -t -o 5 -d Timing_Benchmarks openmp_test128e.fds
$QFDS -t -o 6 -d Timing_Benchmarks openmp_test128f.fds
$QFDS -t -o 7 -d Timing_Benchmarks openmp_test128g.fds
$QFDS -t -o 8 -d Timing_Benchmarks openmp_test128h.fds
