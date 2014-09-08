#!/bin/bash
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh intel64
source ~/.bashrc_fds intel64ib1p4p3

echo Building $target 
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
