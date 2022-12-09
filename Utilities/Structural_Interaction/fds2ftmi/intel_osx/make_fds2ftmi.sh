#!/bin/bash

source $ONEAPI_ROOT/setvars.sh intel64
rm *.o
make -f ../Makefile intel_osx
