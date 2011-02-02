#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source ../SET_MYFDSENV.sh $platform build

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile_openmp $target
