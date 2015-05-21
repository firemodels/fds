#!/bin/bash
source ../setopts.sh $*
LIBDIR=../LIBS/lib_osx_intel_64/
source ../test_libs.sh
KWDIR=../../../Utilities/keyword
SDIR=../../source

source $KWDIR/expand_file.sh $KWDIR $SDIR/smokeview $SDIR/shared/string_util.c
make -f ../Makefile clean
eval make ${SMV_MAKE_OPTS}-f ../Makefile intel_osx_64
source $KWDIR/contract_file.sh $KWDIR $SDIR/smokeview/main.c
