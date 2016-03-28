#!/bin/bash
source ../../scripts/setopts.sh $*
LIBDIR=../../LIBS/lib_osx_intel_64/
source ../../scripts/test_libs.sh

make -f ../Makefile clean
eval make -j 4 ${SMV_MAKE_OPTS}-f ../Makefile intel_osx_64
