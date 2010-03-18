#!/bin/csh -f
set version=$1
set revision=$2
set SVNROOT=~/FDS-SMV

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set LINUXDIR=smv_$version\_linux64

cd $FORBUNDLE

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp readme.html $LINUXDIR/Documentation/release_notes.html
cp $BINDIR/smv5_linux_64 $LINUXDIR/.
cp $BINDIR/smokezip_linux_64 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
