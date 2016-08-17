#!/bin/bash
LIBDIR=../../LIBS/gnu_linux_64/
source ../../scripts/test_libs.sh lua
LUA_SCRIPTING="LUA_SCRIPTING=true"

make -f ../Makefile clean
eval make -j 4 ${LUA_SCRIPTING} -f ../Makefile gnu_linux_64

