#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

if [ "$IFORT_COMPILER" != "" ]; then
source $IFORT_COMPILER/bin/compilervars.sh $platform
fi

echo Building $target
make VPATH="../../../Source:../Source" -f ../makefile $target
