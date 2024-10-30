#!/bin/bash

dir=`pwd`
target=${dir##*/}

# Compile third-party libraries
export FDS_BUILD_TARGET=$target
export SOURCE_INTEL_IFORT=1
source ../Scripts/build_thirdparty_libs.sh "$@"

# Build fds
echo Building $target with Intel MPI and $INTEL_IFORT
make -j4 VPATH="../../Source" -f ../makefile $target
