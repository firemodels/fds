#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source ../Scripts/set_env.sh eth
if [ $? -eq 0 ]; then
  exit
fi

echo Building $target with $MPIDIST
make -j4 MPIFORT="$MPIFORT" VPATH="../../Source" -f ../makefile $target
