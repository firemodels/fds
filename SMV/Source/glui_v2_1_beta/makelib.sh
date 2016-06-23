#!/bin/bash
rm -f *.o
source ../scripts/setopts.sh $*

rm -f *.o
eval make COMPILER=${COMPILER} COMPILER2=${COMPILER2} CFLAGOPT=${CFLAGOPT} SIZE=${SIZE} PLATFORM=\"${PLATFORM}\"
