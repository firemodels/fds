#!/bin/bash

QFDS=../../Utilities/Scripts/qfds.sh  -I 
for i in `seq 1 10`; do
arg=0$i
  if [ $i -gt 9 ]; then
    arg=$i
  fi
  ./makecase64.sh 64 openmp_test64a$arg
  $QFDS openmp_test64a$arg.fds

  ./makecase64.sh 64 openmp_test64d$arg
  $QFDS openmp_test64d$arg.fds

  ./makecase128.sh 128 openmp_test128a$arg
  $QFDS openmp_test128a$arg.fds

  ./makecase128.sh 128 openmp_test128d$arg
  $QFDS openmp_test128d$arg.fds
done
