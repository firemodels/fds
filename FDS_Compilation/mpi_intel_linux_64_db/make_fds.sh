#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh $platform
source ~/.bashrc_fds $platform

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
