#!/bin/bash
platform=intel64
target=intel_linux_64p

source $IFORT_COMPILER/bin/compilervars.sh $platform

rm -f *.o *.mod
make -f ./makefile $target
