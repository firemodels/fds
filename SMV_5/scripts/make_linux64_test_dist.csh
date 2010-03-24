#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set LINUXDIR=smv_test\_$revision\_linux64

cd $FORBUNDLE

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp note.txt $LINUXDIR/Documentation/.
cp $BINDIR/smv5_linux_test_64 $LINUXDIR/.
cp $BINDIR/smokezip_linux_64 $LINUXDIR/.
cp $BINDIR/smokediff_linux_64 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
