#!/bin/csh -f

set directory=$1
set      host=$2
set target=$3

echo Removing smokeview files in $directory on $host
echo
ssh $host \( cd \~/$directory \; make -f ../Makefile $target \)
