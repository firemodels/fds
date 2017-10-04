#!/bin/bash

platform=intel64
dir=`pwd`
target=${dir##*/}

echo Building $target with Intel MPI
make -j4 MPIFORT="$MPIFORT" VPATH="../../Source" -f ../makefile $target
