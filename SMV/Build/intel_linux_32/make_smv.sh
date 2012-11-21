#!/bin/bash
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh ia32
source $IFORT_COMPILER/bin/ifortvars.sh ia32
make -f ../Makefile clean
eval make -j4 ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_32
