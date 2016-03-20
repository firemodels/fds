#!/bin/bash
source $IFORT_COMPILER/bin/iccvars.sh intel64
rm -f *.o
make -f ../Makefile intel_linux_64
