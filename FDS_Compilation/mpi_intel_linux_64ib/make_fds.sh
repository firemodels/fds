#!/bin/bash
platform=ib64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh intel64
source ~/.bashrc_fds $platform

echo Building $target 
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
