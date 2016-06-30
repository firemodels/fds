#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64
rm -f *.o
make -f ../Makefile intel_osx_64
