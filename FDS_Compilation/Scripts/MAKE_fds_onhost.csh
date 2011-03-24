#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3
if ($#argv>3) then
  echo Cleaning FDS in $fdsdir on $host
  ssh $host \( cd \~/$fdsdir \; make -f ../makefile clean \)
  exit 0
endif
echo Building FDS in $fdsdir on $host using the target $target
ssh $host \( cd \~/$fdsdir \; ./make_fds.sh \)
echo Completed
