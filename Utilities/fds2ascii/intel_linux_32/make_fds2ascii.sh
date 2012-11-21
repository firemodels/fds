#!/bin/bash
rm *.o
source $IFORT_COMPILER/bin/ifortvars.sh ia32
make -f ../Makefile intel_linux_32
