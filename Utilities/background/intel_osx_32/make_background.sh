#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh ia32

rm -f *.o
make -f ../Makefile intel_osx_32
