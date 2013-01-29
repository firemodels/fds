#!/bin/bash
# use on systems that do not support -m32
#SIZE=""
SIZE="-m32"
rm -rf *.o
eval make -f ../Makefile SIZE=${SIZE} gcc_linux_32
