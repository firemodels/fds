#!/bin/bash
platform=intel64
target=intel_linux_64

source $IFORT_COMPILER/bin/compilervars.sh $platform

make -f ./makefile $target
