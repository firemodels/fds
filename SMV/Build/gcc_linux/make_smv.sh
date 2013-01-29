#!/bin/bash
# use on systems that do not support -m32
SIZE=""
rm -rf *.o
eval make -f ../Makefile SIZE=${SIZE} gcc_linux
