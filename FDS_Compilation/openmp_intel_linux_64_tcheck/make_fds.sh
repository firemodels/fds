#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source ../SET_FORT.sh $platform
#source /exports/intel/itt/tcheck/bin/32/tcheckvars.csh

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile_openmp $target
