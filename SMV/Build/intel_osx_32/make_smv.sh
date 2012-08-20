#!/bin/bash
source ../setopts.sh $*
make -f ../Makefile clean
eval make ${SMV_MAKE_OPTS} -f ../Makefile intel_osx_32
