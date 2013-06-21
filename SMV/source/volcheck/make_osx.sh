#!/bin/bash
target=intel_osx_64

rm -f *.o *.mod
make -f ./makefile $target
