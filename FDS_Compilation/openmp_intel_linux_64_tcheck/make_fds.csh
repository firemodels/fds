#!/bin/csh -f
set platform=intel64
set dir=`pwd`
#set target=$dir:t
set target=openmp_intel_linux_64_tcheck


source ../Scripts/set_fort.csh $platform
source /exports/intel/itt/tcheck/bin/32/tcheckvars.csh
if ($status == 1) exit

echo Building $target
make VPATH="../../FDS_Source" -f ../makefile_openmp $target
