#!/bin/bash
source ../setopts.sh $*
source ../test_ifort.sh
source $IFORT_COMPILER/bin/compilervars.sh intel64
LIBDIR=../LIBS/lib_linux_intel_64/
source ../test_libs.sh

make -f ../Makefile clean
eval make ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_64_db

