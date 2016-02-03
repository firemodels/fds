#!/bin/bash
rm -f *.o
source ../../setopts.sh $*

rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} ${TARGET} SYSLIBS="-Wl,-E -ldl" CFLAGS="-DLUA_USE_LINUX"
