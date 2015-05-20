#!/bin/bash -f
source $IFORT_COMPILER/bin/compilervars.sh intel64
KWDIR=../../keyword
SDIR=../../../SMV/source

rm -f *.o
$KWDIR/expand_file.sh $KWDIR $SDIR/wind2dfs $SDIR/shared/string_util.c
make -j4 -f ../Makefile intel_linux_64
$KWDIR/contract_file.sh $KWDIR $SDIR/shared/string_util.c
