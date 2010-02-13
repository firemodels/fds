#!/bin/csh -f
set version=$1
set revision=$2
set SVNROOT=~/FDS-SMV

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set LINUXDIR=smv_$version\_$revision\_linux

cd $FORBUNDLE

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp note.txt $LINUXDIR/Documentation/.
cp $BINDIR/smv5_linux_32 $LINUXDIR/.
cp $BINDIR/smokezip_linux_32 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
