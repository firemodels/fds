#!/bin/bash
platform=intel64
dir=`pwd`
target=${dir##*/}

source $IFORT_COMPILER/bin/compilervars.sh $platform
source $IFORT_COMPILER/../impi/5.0.0.028/intel64/bin/mpivars.sh
echo Building $target
make VPATH="../" -f ../makefile $target
#mpiifort -o test_mpi_openmp -openmp ../test_mpi_openmp.f90
