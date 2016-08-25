#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64
source ../../scripts/setopts.sh $*
LIBDIR=../../LIBS/intel_osx_64/
source ../../scripts/test_libs.sh

rm -f *.o
make -f ../Makefile intel_osx_64
