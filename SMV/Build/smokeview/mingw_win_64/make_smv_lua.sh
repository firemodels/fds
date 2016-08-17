# !/bin/bash
# source ../../scripts/setopts.sh $*
LIBDIR=../../LIBS/lib_win_mingw_64/
source ../../scripts/test_libs.sh lua
LUA_SCRIPTING="LUA_SCRIPTING=true"

make -f ../Makefile clean
eval make -j 4 ${SMV_MAKE_OPTS} ${LUA_SCRIPTING} -f ../Makefile mingw_win_64
