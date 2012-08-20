#!/bin/bash -f
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
make -f ../Makefile clean
echo *** $SMV_MAKE_OPTS ***
#make -j4 ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_64
make -j4 -f ../Makefile intel_linux_64
