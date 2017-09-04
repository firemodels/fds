#!/bin/bash
dir=`pwd`
target=${dir##*/}

echo Building $target
# GNU make needs to be used on AIX, usually comes as gmake
gmake -j4 VPATH="../../Source" -f ../makefile $target
