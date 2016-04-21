#!/bin/bash
source ../scripts/setopts.sh $*
source $IFORT_COMPILER/bin/compilervars.sh intel64
rm -f *.o
STDINT="-DHAVE_STDINT_H"
eval make STDINT=${STDINT} COMPILER=${COMPILER} SIZE=${SIZE} libgd.a
