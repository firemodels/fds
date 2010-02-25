#!/bin/csh -f
set platform=intel64
set dir=`pwd`
set target=$dir:t

source /opt/intel/11/bin/ifortvars.csh $platform
source ../Scripts/set_mpi64.csh

echo Building $target
make VPATH="../../FDS_Source" -f ../makefile $target
