#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh $platform
if [[ $MPIDIST != *ib ]] ; then
# define MPIDIST
source ../Scripts/set_mpidist.sh /shared/openmpi_64ib
fi
if [ "$MPIDIST" == "" ]; then
# if MPIDIST was not defined above, abort
  exit
fi

echo Building $target 
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
