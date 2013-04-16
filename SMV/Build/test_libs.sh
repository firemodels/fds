#!/bin/bash -f
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
if [ "`uname`" != "Darwin" ]; then
makelibs libglut.a 
fi
makelibs libjpeg.a 
makelibs libpng.a 
makelibs libz.a 

