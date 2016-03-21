#!/bin/bash
revision=$1
SVNROOT=~/$2

BACKGROUNDDIR=$SVNROOT/SMV/Build/background/intel_linux_64
SMVDIR=$SVNROOT/SMV/Build/smokeview/intel_linux_64
SMZDIR=$SVNROOT/SMV/Build/smokezip/intel_linux_64
SMDDIR=$SVNROOT/SMV/Build/smokediff/intel_linux_64
WINDDIR=$SVNROOT/SMV/Build/wind2fds/intel_linux_64
FORBUNDLE=$SVNROOT/SMV/for_bundle
LINUXDIR=smv_test\_$revision\_linux64
UPDATER=$SVNROOT/Utilities/Scripts/make_updater.sh

cd $SVNROOT/SMV/uploads

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/bin
mkdir -p $LINUXDIR/Documentation
cp -r $FORBUNDLE/textures $LINUXDIR/bin/.
cp $FORBUNDLE/objects.svo $LINUXDIR/bin/.
cp $FORBUNDLE/smokeview.ini $LINUXDIR/bin/.
cp $FORBUNDLE/volrender.ssf $LINUXDIR/bin/.
cp $FORBUNDLE/*.po $LINUXDIR/bin/.
cp $BACKGROUNDDIR/background $LINUXDIR/bin/.
cp $SMVDIR/smokeview_linux_test_64 $LINUXDIR/bin/smokeview
cp $SMZDIR/smokezip_linux_64 $LINUXDIR/bin/smokezip
cp $SMDDIR/smokediff_linux_64 $LINUXDIR/bin/smokediff
cp $WINDDIR/wind2fds_linux_64 $LINUXDIR/bin/wind2fds
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
cd $LINUXDIR
tar cvf ../$LINUXDIR.tar .
cd ..
gzip $LINUXDIR.tar
$UPDATER Linux 64 $revision $LINUXDIR.tar.gz $LINUXDIR.sh FDS/FDS6 test
