#!/bin/bash
fdsrun=../../FDS_Compilation/openmp_intel_linux_64/fds_openmp_intel_linux_64 

qfds.sh -t -o 1 $fdsrun openmp_test64a.fds
qfds.sh -t -o 2 $fdsrun openmp_test64b.fds
qfds.sh -t -o 3 $fdsrun openmp_test64c.fds
qfds.sh -t -o 4 $fdsrun openmp_test64d.fds
qfds.sh -t -o 5 $fdsrun openmp_test64e.fds
qfds.sh -t -o 6 $fdsrun openmp_test64f.fds
qfds.sh -t -o 7 $fdsrun openmp_test64g.fds
qfds.sh -t -o 8 $fdsrun openmp_test64h.fds

qfds.sh -t -o 1 $fdsrun openmp_test128a.fds
qfds.sh -t -o 2 $fdsrun openmp_test128b.fds
qfds.sh -t -o 3 $fdsrun openmp_test128c.fds
qfds.sh -t -o 4 $fdsrun openmp_test128d.fds
qfds.sh -t -o 5 $fdsrun openmp_test128e.fds
qfds.sh -t -o 6 $fdsrun openmp_test128f.fds
qfds.sh -t -o 7 $fdsrun openmp_test128g.fds
qfds.sh -t -o 8 $fdsrun openmp_test128h.fds
