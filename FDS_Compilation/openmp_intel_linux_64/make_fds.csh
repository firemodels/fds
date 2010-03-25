#!/bin/csh -f
set platform=intel64
set dir=`pwd`
set target=openmp_intel_linux_64


source ../Scripts/set_fort.csh $platform
if ($status == 1) exit

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile_openmp $target
