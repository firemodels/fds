#!/bin/bash
dir=`pwd`
target=${dir##*/}

../Scripts/save_fdsinfo.sh
echo Building $target
make -j4 VPATH="../../Source" -f ../makefile $target
