#!/bin/bash

source $$ONEAPI_ROOT/setvars.sh ia32
rm *.o
make -f ../Makefile intel_osx_32
