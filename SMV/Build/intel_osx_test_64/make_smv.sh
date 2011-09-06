#!/bin/bash
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
make -f ../Makefile intel_osx_test_64
