#!/bin/bash
platform=ia32
dir=`pwd`
target=${dir##*/}

source ../SET_FORT.sh $platform
source ~/.bashrc_fds $platform

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
