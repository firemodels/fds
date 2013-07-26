#!/bin/bash -f
source $IFORT_COMPILER/bin/iccvars.sh ia32
make -j4 -f ../Makefile intel_linux_32
