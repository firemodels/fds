#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3
set out=$target.out
if ($#argv>3) then
  echo Cleaning FDS in $fdsdir on $host |& tee $out
  ssh $host \( cd \~/$fdsdir \; make -f ../makefile clean \) |& tee -a $out
  scp $out $host\:$fdsdir/.
  exit 0
endif
echo | tee -a $out
echo Building FDS in $fdsdir on $host using the target $target |& tee -a $out
ssh $host \( cd \~/$fdsdir \; make VPATH="../../../FDS_Source" -f ../makefile $target \) |& tee -a $out
scp $out $host\:$fdsdir/.
