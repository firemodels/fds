#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3
if ($#argv>3) then
  echo Cleaning FDS in $fdsdir
  ssh $host \( cd \~/$fdsdir \; make -f ../makefile clean \)
  exit 0
endif
echo Building FDS using:
echo    target: $target
echo directory: $fdsdir
echo      host: $host
echo
ssh $host \( cd \~/$fdsdir \; make VPATH="../../../FDS_Source" -f ../makefile $target \)
