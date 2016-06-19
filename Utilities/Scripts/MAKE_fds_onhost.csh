#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3
set out=~/$fdsdir/$target.out
if ($#argv>3) then
  date |& tee $out
  echo Cleaning FDS in $fdsdir on $host |& tee -a $out
  ssh $host \( cd \~/$fdsdir \; make -f ../makefile clean \) |& tee -a $out
  exit 0
endif
echo |& tee -a $out
date |& tee -a $out
echo Building FDS in $fdsdir on $host using the target $target |& tee -a $out
ssh $host svn info $fdsdir/../../../FDS_Source/ | grep Rev | grep Last |& tee -a $out
ssh $host \( cd \~/$fdsdir \; make VPATH="../../../FDS_Source" -f ../makefile $target \) |& tee -a $out
echo Completed |& tee -a $out
