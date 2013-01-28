#!/bin/bash
rm -f *.o
#eval make -j4 ${SMV_MAKE_OPTS} 
COMPILER=gcc
#COMPILER=icc
#SIZE="-m32"
SIZE="-m64"
eval make COMPILER=${COMPILER} SIZE=${SIZE} libgd.a
