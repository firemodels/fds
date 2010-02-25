#!/bin/csh -f
set platform=intel64
set dir=`pwd`
set target=$dir:t

source /opt/intel/11/bin/ifortvars.csh $platform

echo Building $target
make VPATH="../../FDS_Source" -f ../makefile $target
