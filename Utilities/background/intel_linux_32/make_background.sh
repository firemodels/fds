#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh ia32
KWDIR=../../keyword
SDIR=../../../SMV/source

rm -f *.o
$KWDIR/expand_file.sh $KWDIR $SDIR/background $SDIR/shared/string_util.c
make -j4 -f ../Makefile intel_linux_32
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
