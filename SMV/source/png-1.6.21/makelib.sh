#!/bin/bash
rm -f *.o
source ../scripts/setopts.sh $*
source $IFORT_COMPILER/bin/compilervars.sh intel64

rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} libpng.a
