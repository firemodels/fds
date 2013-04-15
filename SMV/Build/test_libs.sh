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
makelibs libglut.a 
makelibs libjpeg.a 
makelibs libpng.a 
makelibs libz.a 

