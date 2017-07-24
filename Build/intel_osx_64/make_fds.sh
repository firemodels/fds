#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

if [ "$IFORT_COMPILER" != "" ]; then
  source $IFORT_COMPILER/bin/compilervars.sh $platform
fi
VERSION=`ifort -v 2>&1 | awk '{print $3}' | awk -F'.' '{print $1}'`


echo Building $target
if [ "$VERSION" == "16" ]; then
  make -j4 FOPENMPFLAGS="-openmp -openmp-link static" VPATH="../../Source" -f ../makefile $target
else
  make -j4 VPATH="../../Source" -f ../makefile $target
fi
