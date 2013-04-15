#!/bin/bash -f
source ../setopts.sh $*
source $IFORT_COMPILER/bin/iccvars.sh intel64
source $IFORT_COMPILER/bin/ifortvars.sh intel64
LIBDIR=../LIBS/lib_linux_intel_64/

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
eval make -j4 ${SMV_MAKE_OPTS} -f ../Makefile intel_linux_64

