#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh ia32
KWDIR=../../keyword
SDIR=../../../SMV/source

$KWDIR/expand_file.sh $KWDIR $SDIR/background $SDIR/shared/string_util.c
rm -f *.o
make -j4 -f ../Makefile intel_osx_32
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
