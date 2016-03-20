#!/bin/bash
source $IFORT_COMPILER/bin/iccvars.sh intel64
rm -f *.o
make -j4 -f ../Makefile intel_linux_64_db
