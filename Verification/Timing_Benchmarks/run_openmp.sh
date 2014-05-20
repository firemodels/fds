#!/bin/bash
fdsrun=../../FDS_Compilation/openmp_intel_linux_64/fds_openmp_intel_linux_64 
./makecase.sh test1
./makecase.sh test2
./makecase.sh test3
./makecase.sh test4
./makecase.sh test5
./makecase.sh test6
./makecase.sh test7
./makecase.sh test8

qfds.sh -o 1 $fdsrun test1.fds
qfds.sh -o 2 $fdsrun test2.fds
qfds.sh -o 3 $fdsrun test3.fds
qfds.sh -o 4 $fdsrun test4.fds
qfds.sh -o 5 $fdsrun test5.fds
qfds.sh -o 6 $fdsrun test6.fds
qfds.sh -o 7 $fdsrun test7.fds
qfds.sh -o 8 $fdsrun test8.fds
