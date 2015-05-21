#!/bin/bash
KWDIR=../../keyword
SDIR=../../../SMV/source

rm -f *.o
$KWDIR/expand_file.sh $KWDIR $SDIR/wind2dfs $SDIR/shared/string_util.c
make -f ../Makefile intel_osx_64
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
