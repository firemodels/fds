#!/bin/bash
QFDS=../../Utilities/Scripts/qfds.sh

file=openmp64_test0
./makecase.sh 64 $file
$QFDS -t -o 1 -d . $file.fds
for index in {1..20}
do
  file=openmp64_test$index
  ./makecase.sh 64 $file
  $QFDS -t -o 4 -d . $file.fds
done

file=openmp128_test0
./makecase.sh 128 $file
$QFDS -t -o 1 -d . $file.fds
for index in {1..20}
do
  file=openmp128_test$index
  ./makecase.sh 128 $file
  $QFDS -t -o 4 -d . $file.fds
done
