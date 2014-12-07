#!/bin/bash
QFDS=../../Utilities/Scripts/qfds.sh

for index in {1..50}
do
  file=openmp_test$index
  ./makecase.sh 64 $file
  $QFDS -t -o 4 -d . $file.fds
done
