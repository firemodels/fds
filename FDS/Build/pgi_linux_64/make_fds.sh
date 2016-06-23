#!/bin/bash
dir=`pwd`
target=${dir##*/}


echo Building $target
make VPATH="../../Source" -f ../makefile $target
