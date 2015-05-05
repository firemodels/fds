#!/bin/bash
source ../setopts.sh $*
#source ../test_ifort.sh
LIBDIR=../LIBS/lib_osx_intel_64/
source ../test_libs.sh

eval make ${SMV_MAKE_OPTS}-f ../Makefile intel_osx_64
