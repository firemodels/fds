#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64
rm *.o
make -j4 -f ../Makefile intel_linux_64
