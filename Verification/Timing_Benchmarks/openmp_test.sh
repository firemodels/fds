#!/bin/bash
QFDS=../../Utilities/Scripts/qfds.sh

for index in {1..10}
do
  file=openmp64_test1_$index
  ./makecase.sh 64 $file
  $QFDS -t -o 1 -d . $file.fds
done

for index in {1..10}
do
  file=openmp64_test4_$index
  ./makecase.sh 64 $file
  $QFDS -t -o 4 -d . $file.fds
done
