#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV

set BACKGROUNDDIR=$SVNROOT/Utilities/background/intel_linux_32
set SMVDIR=$SVNROOT/SMV/Build/intel_linux_32
set SMZDIR=$SVNROOT/Utilities/smokezip/intel_linux_32
set SMDDIR=$SVNROOT/Utilities/smokediff/intel_linux_32
set WINDDIR=$SVNROOT/Utilities/wind2fds/intel_linux_32
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set LINUXDIR=smv_test\_$revision\_linux32
set UPDATER=$SVNROOT/Utilities/Scripts/make_updater.sh

cd $FORBUNDLE/uploads

rm -rf $LINUXDIR
mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/bin
mkdir -p $LINUXDIR/Documentation
cp $FORBUNDLE/volrender.ssf $LINUXDIR/bin/.
cp $FORBUNDLE/objects.svo $LINUXDIR/bin/.
cp $FORBUNDLE/*.po $LINUXDIR/bin/.
cp $BACKGROUNDDIR/background $LINUXDIR/bin/.
cp $SMVDIR/smokeview_linux_test_32 $LINUXDIR/bin/smokeview
cp $SMZDIR/smokezip_linux_32 $LINUXDIR/bin/smokezip
cp $SMDDIR/smokediff_linux_32 $LINUXDIR/bin/smokediff
cp $WINDDIR/wind2fds_linux_32 $LINUXDIR/bin/wind2fds
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
cd $LINUXDIR
tar cvf ../$LINUXDIR.tar .
cd ..
gzip $LINUXDIR.tar
$UPDATER Linux 32 $revision $LINUXDIR.tar.gz $LINUXDIR.sh FDS/FDS6 test
