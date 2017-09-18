#!/bin/bash

dir=`pwd`
target=${dir##*/}

echo Building $target with $MPIDIST
`ifort -v`

make -j4 VPATH="../../Source" -f ../makefile $target
