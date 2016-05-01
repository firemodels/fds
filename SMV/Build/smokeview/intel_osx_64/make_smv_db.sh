#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64

make -f ../Makefile clean
make -j 4 -f ../Makefile intel_osx_64_db
