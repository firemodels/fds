#!/bin/bash
rm *.o
source $IFORT_COMPILER/bin/ifortvars.sh intel64
make -f ../Makefile intel_linux_64
