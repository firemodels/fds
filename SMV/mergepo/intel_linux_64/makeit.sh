#!/bin/bash
source $IFORT_COMPILER/bin/iccvars.sh intel64
rm *.o
make -f ../Makefile intel_linux_64
