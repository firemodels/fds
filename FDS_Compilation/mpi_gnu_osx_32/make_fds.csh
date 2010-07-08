#!/bin/csh -f
set platform=intel32
set dir=`pwd`
set target=$dir:t


source ../Scripts/set_gnufort_osx.csh $platform
if ($status == 1) exit

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
