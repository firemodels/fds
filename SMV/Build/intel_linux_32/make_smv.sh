#!/bin/bash
source ../clean.sh $*
source $IFORT_COMPILER/bin/iccvars.sh ia32
source $IFORT_COMPILER/bin/ifortvars.sh ia32
make -f ../Makefile clean
make -j4 SMV_TESTFLAG="$SMV_TESTFLAG" SMV_TESTSTRING="$SMV_TESTSTRING" -f ../Makefile intel_linux_32
