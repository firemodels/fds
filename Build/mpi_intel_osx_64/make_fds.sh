#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}
source $IFORT_COMPILER/bin/compilervars.sh $platform
source ../Scripts/set_mpidist.sh eth $MPIDIST_ETH

echo Building $target
make -j4 VPATH="../../Source" -f ../makefile $target
../../Utilities/Scripts/md5hash.sh fds_mpi_intel_osx_64
