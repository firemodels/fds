#!/bin/bash

# build using 64 bit intel compilers
OPTS="-i -6"
source $IFORT_COMPILER/bin/compilervars.sh intel64

LIBDIR=`pwd`
SRCDIR=$LIBDIR/../../../source
cd $SRCDIR
SRCDIR=`pwd`

# GD
cd $SRCDIR/gd-2.0.15
./makelib.sh $OPTS
cp libgd.a $LIBDIR/.

# GLUI
cd $SRCDIR/glui_v2_1_beta
./makelib.sh $OPTS
cp libglui.a $LIBDIR/.

# GLUT
cd $SRCDIR/glut-3.7.6
./makelib.sh $OPTS
cp libglut.a $LIBDIR/.

# JPEG
cd $SRCDIR/jpeg-6b
./makelib.sh $OPTS
cp libjpeg.a $LIBDIR/.

# PNG
cd $SRCDIR/png125
./makelib.sh $OPTS
cp libpng.a $LIBDIR/.

# ZLIB
cd $SRCDIR/zlib114
./makelib.sh $OPTS
cp libz.a $LIBDIR/.
