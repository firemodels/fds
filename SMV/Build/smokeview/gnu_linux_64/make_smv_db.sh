#!/bin/bash
LIBDIR=../../LIBS/gnu_linux_64/
source ../../scripts/test_libs.sh

make -f ../Makefile clean
eval make -j 4 -f ../Makefile gnu_linux_64_db

