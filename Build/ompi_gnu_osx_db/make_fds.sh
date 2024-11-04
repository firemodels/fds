#!/bin/bash

dir=`pwd`
target=${dir##*/}

# Compile third-party libraries
export FDS_BUILD_TARGET=$target
source ../Scripts/build_thirdparty_libs.sh "$@"

# Build fds
echo Building $target
make -j4 VPATH="../../Source" -f ../makefile $target
