#!/bin/bash

# Determine if integers in last row (col column) of files
# file1 and file2 are with tol of each other.
# Return 1 if yes, 0 if no

file1=$1
file2=$2
col=$3
tol=$4

if ! [ -e $file1 ]; then
  echo "*** Warning: file $file1 does not exist"
  echo "command: $0 $*"
  exit
fi
if ! [ -e $file2 ]; then
  echo "*** Warning: file $file2 does not exist"
  echo "command: $0 $*"
  exit
fi

 # Absolute value
 abs() {
   [ $1 -lt 0 ] && echo $((-$1)) || echo $1
 }

num1=`tail -n 1 $file1 | awk -v var="$col" -F',' '{print $var}'`
num2=`tail -n 1 $file2 | awk -v var="$col" -F',' '{print $var}'`
diff=`expr $num1 - $num2`
diff=`abs $diff`
if [ $diff -gt $tol ]; then
  echo "*** Warning: column $col of file $file2 is out of tolerance"
  echo "             |diff|=$diff, tolerance=$tol"
fi
