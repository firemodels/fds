#!/bin/csh -f
set mssg="64 bit MPI Linux FDS"
set target=mpi_intel_linux_64
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
make VPATH="../../../FDS_Source" -f ../makefile $target |& tee -a $out
date |& tee -a $out
echo Complete |& tee -a $out
