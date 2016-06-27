#!/bin/bash
OPTS="-g -6"
LIBDIR=`pwd`
SRCDIR=$LIBDIR/../../../Source
rm *.a
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
#cd $SRCDIR/glut-3.7.6
#./makelib.sh $OPTS
#cp libglut.a $LIBDIR/.

# JPEG
cd $SRCDIR/jpeg-9b
./makelib.sh $OPTS
cp libjpeg.a $LIBDIR/.

# PNG
cd $SRCDIR/png-1.6.21
./makelib.sh $OPTS
cp libpng.a $LIBDIR/.

# ZLIB
cd $SRCDIR/zlib128
./makelib.sh $OPTS
cp libz.a $LIBDIR/.
