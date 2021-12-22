#!/bin/bash
ARG=$1

source ../Scripts/set_intel_compiler.sh $ARG

dir=`pwd`
target=${dir##*/}

echo Building $target with Intel MPI and $INTEL_IFORT
make VPATH="../../Source" -f ../makefile $target
