#!/bin/bash
source ../scripts/setopts.sh $*
rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} libgd.a
