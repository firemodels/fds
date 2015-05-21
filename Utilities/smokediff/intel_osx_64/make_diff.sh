#!/bin/bash
rm -f *.o
KWDIR=../../keyword
SDIR=../../../SMV/source

$KWDIR/expand_file.sh $KWDIR $SDIR/smokediff $SDIR/shared/string_util.c
make -j4 -f ../Makefile intel_osx_64
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
