#!/bin/bash

source ~/.bashrc_setfort ia32
rm *.o
make -f ../Makefile intel_osx_32
