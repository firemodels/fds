#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3

echo Building FDS using:
echo    target: $target
echo directory: $fdsdir
echo      host: $host
echo
ssh $host \( cd \~/$fdsdir \; make VPATH="../../../FDS_Source" -f ../makefile $target \)
