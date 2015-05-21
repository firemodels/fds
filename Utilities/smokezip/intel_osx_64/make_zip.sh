#!/bin/bash
KWDIR=../../keyword
SDIR=../../../SMV/source

$KWDIR/expand_file.sh $KWDIR $SDIR/smokezip $SDIR/shared/string_util.c
rm -f *.o
make -f ../Makefile intel_osx_64
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
