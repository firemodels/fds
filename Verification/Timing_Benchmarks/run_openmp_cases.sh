#!/bin/bash
fdsrun=../../FDS_Compilation/openmp_intel_linux_64/fds_openmp_intel_linux_64 

qfds.sh -t -o 1 $fdsrun test64a.fds
qfds.sh -t -o 2 $fdsrun test64b.fds
qfds.sh -t -o 3 $fdsrun test64c.fds
qfds.sh -t -o 4 $fdsrun test64d.fds
qfds.sh -t -o 5 $fdsrun test64e.fds
qfds.sh -t -o 6 $fdsrun test64f.fds
qfds.sh -t -o 7 $fdsrun test64g.fds
qfds.sh -t -o 8 $fdsrun test64h.fds

qfds.sh -t -o 1 $fdsrun test128a.fds
qfds.sh -t -o 2 $fdsrun test128b.fds
qfds.sh -t -o 3 $fdsrun test128c.fds
qfds.sh -t -o 4 $fdsrun test128d.fds
qfds.sh -t -o 5 $fdsrun test128e.fds
qfds.sh -t -o 6 $fdsrun test128f.fds
qfds.sh -t -o 7 $fdsrun test128g.fds
qfds.sh -t -o 8 $fdsrun test128h.fds
