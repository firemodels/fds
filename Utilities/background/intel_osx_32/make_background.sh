#!/bin/bash -f
source $IFORT_COMPILER/bin/ifortvars.sh ia32
make -j4 -f ../Makefile intel_osx_32
