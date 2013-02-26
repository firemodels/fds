#!/bin/bash -f
rm *.o
make -j4 -f ../Makefile gcc_linux_64
