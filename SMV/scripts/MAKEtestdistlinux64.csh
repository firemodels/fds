#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV

set SMVDIR=$SVNROOT/SMV/Build/INTEL_LINUX_TEST_64
set SMZDIR=$SVNROOT/Utilities/smokezip/INTEL_LINUX_64
set SMDDIR=$SVNROOT/Utilities/smokediff/INTEL_LINUX_64
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set LINUXDIR=smv_test\_$revision\_linux64

cd $FORBUNDLE

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp note.txt $LINUXDIR/Documentation/.
cp $FORBUNDLE/objects.svo $LINUXDIR/.
cp $FORBUNDLE/*.po $LINUXDIR/.
cp $SMVDIR/smokeview_linux_test_64 $LINUXDIR/.
cp $SMZDIR/smokezip_linux_64 $LINUXDIR/.
cp $SMDDIR/smokediff_linux_64 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
