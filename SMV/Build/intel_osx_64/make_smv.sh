#!/bin/bash
source ../setopts.sh $*

LIBDIR=../LIBS/lib_osx_intel_64/

makelibs()
{
  lib=$1
  if [ ! -e $LIBDIR/$lib ] ; then
    CURDIR=`pwd`
    cd $LIBDIR
    ./makelibs.sh
    cd $CURDIR
  fi
}

makelibs libgd.a
makelibs libglui.a
makelibs libglut.a
makelibs libjpeg.a
makelibs libpng.a
makelibs libz.a

make -f ../Makefile clean
eval make -j4 ${SMV_MAKE_OPTS}-f ../Makefile intel_osx_64
