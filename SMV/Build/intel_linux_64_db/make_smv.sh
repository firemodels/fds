#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64
make -f ../Makefile clean
make -j4 -f ../Makefile intel_linux_64_db
