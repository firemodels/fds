#!/bin/bash
rm -f *.o

make "FORTLIBDIR=$IFORT_COMPILER_LIB" -f ../Makefile intel_osx_64
