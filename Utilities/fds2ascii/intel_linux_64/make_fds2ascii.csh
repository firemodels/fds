#!/bin/csh -f

source $IFORT_COMPILER11/bin/ifortvars.csh intel64

make -f ../Makefile intel_linux_64
