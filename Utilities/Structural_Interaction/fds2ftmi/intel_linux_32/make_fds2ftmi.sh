#!/bin/bash
rm *.o
source $ONEAPI_ROOT/setvars.sh ia32
make -f ../Makefile intel_linux_32
