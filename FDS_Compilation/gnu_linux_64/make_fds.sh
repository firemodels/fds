#!/bin/bash
dir=`pwd`
target=${dir##*/}

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
