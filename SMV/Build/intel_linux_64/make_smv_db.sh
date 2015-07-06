#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64

KWDIR=../../../Utilities/keyword
SDIR=../../source

make -f ../Makefile clean
make -f ../Makefile intel_linux_64_db
