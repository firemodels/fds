#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV

set BACKGROUNDDIR=$SVNROOT/Utilities/background/intel_linux_32
set SMVDIR=$SVNROOT/SMV/Build/intel_linux_64
set SMZDIR=$SVNROOT/Utilities/smokezip/intel_linux_64
set SMDDIR=$SVNROOT/Utilities/smokediff/intel_linux_64
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set LINUXDIR=smv_test\_$revision\_linux64

cd $FORBUNDLE/uploads

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp $FORBUNDLE/objects.svo $LINUXDIR/.
cp $FORBUNDLE/*.po $LINUXDIR/.
cp $BACKGROUNDDIR/background $LINUXDIR/.
cp $SMVDIR/smokeview_linux_test_64 $LINUXDIR/.
cp $SMZDIR/smokezip_linux_64 $LINUXDIR/.
cp $SMDDIR/smokediff_linux_64 $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
