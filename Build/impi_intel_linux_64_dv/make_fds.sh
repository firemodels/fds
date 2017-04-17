#!/bin/bash

platform=intel64
dir=`pwd`
target=${dir##*/}

if [ "$IFORT_COMPILER" != "" ]; then
source $IFORT_COMPILER/bin/compilervars.sh $platform
fi

echo Building $target with Intel MPI
make -j4 MPIFORT="$MPIFORT" VPATH="../../Source" -f ../makefile $target
