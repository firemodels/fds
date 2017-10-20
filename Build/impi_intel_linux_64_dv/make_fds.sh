#!/bin/bash

dir=`pwd`
target=${dir##*/}

echo Building $target with Intel MPI
make -j4 VPATH="../../Source" -f ../makefile $target
