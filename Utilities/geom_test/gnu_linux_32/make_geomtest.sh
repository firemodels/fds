#!/bin/bash
dir=`pwd`
target=${dir##*/}

echo Building $target
make -j4 VPATH="../../../Source:../Source" -f ../makefile $target
