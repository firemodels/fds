#!/bin/bash
source ../../scripts/setopts.sh $*
source ../../scripts/test_ifort.sh
source $IFORT_COMPILER/bin/compilervars.sh intel64
LIBDIR=../../LIBS/intel_linux_64/
source ../../scripts/test_libs.sh

rm -f *.o
make -f ../Makefile intel_linux_64
