#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh $platform
source $I_MPI_ROOT/bin/mpivars.sh
source ~/.bashrc_fds intel64intel

echo Building $target
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
