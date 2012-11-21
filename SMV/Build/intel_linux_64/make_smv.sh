#!/bin/bash -f
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
make -f ../Makefile clean
eval make -j4 ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_64
