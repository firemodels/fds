#!/bin/bash
source ../setopts.sh $*
make -f ../Makefile clean
eval make -j4 ${SMV_MAKE_OPTS}-f ../Makefile intel_osx_64
