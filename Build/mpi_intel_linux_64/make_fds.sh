#!/bin/bash
dir=`pwd`
target=${dir##*/}

echo Building $target
make -j4 MPIFORT="$MPIFORT" VPATH="../../Source" -f ../makefile $target
