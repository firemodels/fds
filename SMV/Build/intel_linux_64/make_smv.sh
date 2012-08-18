#!/bin/bash -f
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
make -f ../Makefile clean
make -j4 SMV_TESTFLAG="$SMV_TESTFLAG" SMV_TESTSTRING="$SMV_TESTSTRING" SMV_PROFILEFLAG="$SMV_PROFILEFLAG" SMV_PROFILESTRING="$SMV_PROFILESTRING" -f ../Makefile intel_linux_64
