#!/bin/csh -f
set mssg="64 bit MPI OSX FDS"
set target=mpi_intel_osx_64
set objdir=.
if($#argv>0)then
  set objdir=$1
endif
cd $objdir

if ($#argv > 1) then
echo Cleaning $mssg
make  -f ../makefile clean
exit 0
endif
echo Building $mssg
make VPATH="../../FDS_Source" -f ../makefile $target
echo Complete

