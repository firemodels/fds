#!/bin/csh -f
set platform=intel64
set dir=`pwd`
set target=$dir:t


source ../SET_MYFDSENV.sh ib64
if ($status == 1) exit

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
