#!/bin/bash
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh ia32
source $IFORT_COMPILER/bin/ifortvars.sh ia32
LIBDIR=../LIBS/lib_linux_intel_32/

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
eval make -j4 ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_32
