#!/bin/bash
rm -f *.o
source ../scripts/setopts.sh $*

rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} libjpeg.a
