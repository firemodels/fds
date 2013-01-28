#!/bin/bash
rm -f *.o

#COMPILER=gcc
COMPILER=icc

#SIZE="-m32"
SIZE="-m64"

rm -f *.o
eval make COMPILER=${COMPILER} SIZE=${SIZE} libjpeg.a
