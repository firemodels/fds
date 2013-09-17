#!/bin/bash
rm *.o
source $IFORT_COMPILER/bin/compilervars.sh intel64
make -f ../Makefile intel_linux_64
