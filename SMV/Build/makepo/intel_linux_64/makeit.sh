#!/bin/bash
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
rm *.o
make -j4 -f ../Makefile intel_linux_64
