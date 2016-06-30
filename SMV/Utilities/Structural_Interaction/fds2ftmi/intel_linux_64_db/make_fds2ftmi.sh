#!/bin/bash
rm -f *.o
source $IFORT_COMPILER/bin/compilervars.sh intel64
make -f ../Makefile intel_linux_64_db
