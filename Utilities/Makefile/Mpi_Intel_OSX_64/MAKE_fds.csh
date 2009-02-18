#!/bin/csh -f
 
cd $1
if ($#argv > 1) then
echo Cleaning 64 bit paraallel OSX FDS
make  -f ../makefile clean
exit 0
endif
make VPATH="../../../FDS_Source" -f ../makefile mpi_intel_osx_64
