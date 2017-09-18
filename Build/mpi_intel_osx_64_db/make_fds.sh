#!/bin/bash

dir=`pwd`
target=${dir##*/}

if [ "$arg" != "-f" ]; then
  source ../Scripts/set_mpidist.sh eth
  if [ $mpi_error -eq 1 ]; then
    exit
  fi
  if [ $setup_fortran -eq 1 ]; then
    source $IFORT_COMPILER/bin/compilervars.sh intel64
  fi
fi
echo Building $target with $MPIDIST
`ifort -v`

make -j4 VPATH="../../Source" -f ../makefile $target
