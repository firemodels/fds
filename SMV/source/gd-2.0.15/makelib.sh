#!/bin/bash
source ../setopts.sh $*
rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} libgd.a
