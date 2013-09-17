#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh ia32
rm *.o
make -j4 -f ../Makefile intel_linux_32
