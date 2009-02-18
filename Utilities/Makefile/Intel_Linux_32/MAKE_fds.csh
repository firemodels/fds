#!/bin/csh -f
set mssg="32 bit Linux FDS"
set target=intel_linux_32

set out=$target.out

cd $1
if ($#argv > 1) then
date | tee $out
echo Cleaning $mssg |& tee -a $out
make  -f ../makefile clean |& tee -a $out
exit 0
endif
echo "" |& tee -a $out
date |& tee -a $out
echo Building $mssg |& tee -a $out
svn info ../../../FDS_Source/ | grep Rev | grep Last |& tee -a $out
make VPATH="../../../FDS_Source" -f ../makefile $target |& tee -a $out
date |& tee -a $out
echo Complete |& tee -a $out
