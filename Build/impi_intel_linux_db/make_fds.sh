#!/bin/bash
ARG=$1

source ../Scripts/set_intel_compiler.sh $ARG

dir=`pwd`
target=${dir##*/}

# build hypre
source ../Scripts/HYPRE/build_hypre.sh confmake_impi_intel_linux.sh

## build sundials
#source ../Scripts/SUNDIALS/build_sundials.sh

# build fds
echo Building $target with Intel MPI and $INTEL_IFORT
make VPATH="../../Source" -f ../makefile $target
