#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/ifortvars.sh $platform

echo Building $target
make VPATH="../../FDS_Source" -f ../makefile $target
