#!/bin/bash

dir=`pwd`
target=${dir##*/}

echo Building $target with Intel MPI
make VPATH="../../Source" -f ../makefile $target
