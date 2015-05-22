#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}
KWDIR=../../Utilities/keyword
SDIR=../../FDS_Source

source $IFORT_COMPILER/bin/compilervars.sh $platform

echo Building $target
source $KWDIR/expand_file.sh $KWDIR $SDIR $SDIR/main.f90
make -j4 VPATH="../../FDS_Source" -f ../makefile $target
source $KWDIR/contract_file.sh $KWDIR $SDIR/main.f90
