#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh ia32
rm *.o
make -f ../Makefile intel_osx_32
