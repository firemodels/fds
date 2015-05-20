#!/bin/bash
source $IFORT_COMPILER/bin/compilervars.sh intel64
KWDIR=../../keyword
SDIR=../../../SMV/source

rm -f *.o
$KWDIR/expand_file.sh $KWDIR $SDIR/smokediff $SDIR/shared/string_util.c
make -f ../Makefile intel_linux_64
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
