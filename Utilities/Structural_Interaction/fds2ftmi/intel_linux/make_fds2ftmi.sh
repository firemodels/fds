#!/bin/bash
rm *.o
source $ONEAPI_ROOT/setvars.sh intel64
make -f ../Makefile intel_linux
