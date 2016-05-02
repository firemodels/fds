#!/bin/bash
DIR="$( cd "$( dirname "${BASH_SOURCE[0]}" )" && pwd )"
cd $DIR
rm *.a

LUA=$1

# build using 64 bit intel compilers
OPTS="-i -6"
source $IFORT_COMPILER/bin/compilervars.sh intel64

LIBDIR=`pwd`
SRCDIR=$LIBDIR/../../../source
cd $SRCDIR
SRCDIR=`pwd`

# GLUT
cd $SRCDIR/glut-3.7.6
./makelib.sh $OPTS
cp libglut.a $LIBDIR/.

# GLUI
cd $SRCDIR/glui_v2_1_beta
./makelib.sh $OPTS
cp libglui.a $LIBDIR/.

# ZLIB
cd $SRCDIR/zlib128
./makelib.sh $OPTS
cp libz.a $LIBDIR/.

# JPEG
cd $SRCDIR/jpeg-9b
./makelib.sh $OPTS
cp libjpeg.a $LIBDIR/.

# PNG
cd $SRCDIR/png-1.6.21
./makelib.sh $OPTS
cp libpng.a $LIBDIR/.

# GD
cd $SRCDIR/gd-2.1.1
./makelib.sh $OPTS
cp libgd.a $LIBDIR/.

if [ "$LUA" == "lua" ]; then

# Lua # Lua interpreter
cd $SRCDIR/lua-5.3.1/src
export TARGET=liblua.a
./makelib.sh $OPTS
cp liblua.a $LIBDIR/.

# LPEG # Lua parsing libarary to parse SSF files
cd $SRCDIR/lpeg-1.0.0
export TARGET=linux
./makelib.sh $OPTS
cp lpeg.so $LIBDIR/.
fi
