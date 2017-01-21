#!/bin/bash

platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh $platform
source ../Scripts/set_mpidist.sh eth $MPIDIST_ETH

if [ "$MPIDIST" == "" ]; then
  exit
fi

echo Building $target with $MPIDIST
make -j4 MPIFORT="$MPIFORT" VPATH="../../Source" -f ../makefile $target
../../Utilities/Scripts/md5hash.sh fds_mpi_intel_linux_64
