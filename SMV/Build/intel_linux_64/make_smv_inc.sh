#!/bin/bash
source ../setopts.sh $*
source ../test_ifort.sh
source $IFORT_COMPILER/bin/compilervars.sh intel64
LIBDIR=../LIBS/lib_linux_intel_64/
source ../test_libs.sh

source $KWDIR/expand_file.sh $KWDIR $SDIR/smokeview $SDIR/shared/string_util.c
eval make ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_64
source $KWDIR/contract_file.sh $KWDIR $SDIR/smokeview/main.c
