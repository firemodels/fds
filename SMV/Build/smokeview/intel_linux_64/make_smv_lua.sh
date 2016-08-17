#!/bin/bash
source ../../scripts/setopts.sh $*
source ../../scripts/test_ifort.sh
source $IFORT_COMPILER/bin/compilervars.sh intel64
LIBDIR=../../LIBS/intel_linux_64/
source ../../scripts/test_libs.sh lua
LUA_SCRIPTING="LUA_SCRIPTING=true"

make -f ../Makefile clean
eval make -j 4 ${SMV_MAKE_OPTS} ${LUA_SCRIPTING} -f ../Makefile intel_linux_64

