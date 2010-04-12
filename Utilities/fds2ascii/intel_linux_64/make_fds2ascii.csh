#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh intel64
rm *.o
make -f ../Makefile intel_linux_64
