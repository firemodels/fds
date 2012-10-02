#!/bin/bash
platform=ia32
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/ifortvars.sh $platform

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
