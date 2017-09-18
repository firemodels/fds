#!/bin/bash

dir=`pwd`
target=${dir##*/}

echo Building $target with $MPIDIST
echo $IFORT_COMPILER
`ifort -v`

make -j4 VPATH="../../Source" -f ../makefile $target
