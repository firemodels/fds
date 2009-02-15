#!/bin/csh -f
 
cd $1

make VPATH="../../../FDS_Source" -f ../makefile mpi_intel_osx_32
