#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV

set BINDIR=$SVNROOT/SMV/bin
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set LINUXDIR=smv_test\_$revision\_linux

cd $FORBUNDLE

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp note.txt $LINUXDIR/Documentation/.
cp $FORBUNDLE/objects.svo $LINUXDIR/.
cp $BINDIR/smv5_linux_test_32 $LINUXDIR/.
cp $BINDIR/smokezip_linux_32 $LINUXDIR/.
cp $BINDIR/smokediff_linux_32 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
