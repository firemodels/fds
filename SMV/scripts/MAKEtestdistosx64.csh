#!/bin/csh -f
set revision=$1
set REMOTESVNROOT=FDS-SMV
set OSXHOST=$2
set SVNROOT=~/$3

set BACKGROUNDDIR=$REMOTESVNROOT/Utilities/background/intel_osx_32
set SMVDIR=$REMOTESVNROOT/SMV/Build/intel_osx_64
set SMZDIR=$REMOTESVNROOT/Utilities/smokezip/intel_osx_64
set SMDDIR=$REMOTESVNROOT/Utilities/smokediff/intel_osx_64
set WINDDIR=$REMOTESVNROOT/Utilities/wind2fds/intel_osx_64
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set OSXDIR=smv_test\_$revision\_osx64
set UPDATER=$SVNROOT/Utilities/Scripts/make_updater.sh

cd $SVNROOT/SMV/uploads

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/bin
mkdir -p $OSXDIR/Documentation

cp $FORBUNDLE/objects.svo $OSXDIR/bin/.
cp $FORBUNDLE/smokeview.ini $OSXDIR/bin/.
cp -r $FORBUNDLE/textures $OSXDIR/bin/.
cp $FORBUNDLE/*.po $OSXDIR/bin/.
cp $FORBUNDLE/volrender.ssf $OSXDIR/bin/.
scp $OSXHOST\:$BACKGROUNDDIR/background $OSXDIR/bin/.
scp $OSXHOST\:$SMVDIR/smokeview_osx_test_64 $OSXDIR/bin/smokeview
scp $OSXHOST\:$SMZDIR/smokezip_osx_64 $OSXDIR/bin/smokezip
scp $OSXHOST\:$SMDDIR/smokediff_osx_64 $OSXDIR/bin/smokediff
scp $OSXHOST\:$WINDDIR/wind2fds_osx_64 $OSXDIR/bin/wind2fds
rm -f $OSXDIR.tar $OSXDIR.tar.gz
cd $OSXDIR
tar cvf ../$OSXDIR.tar .
cd ..
gzip $OSXDIR.tar
$UPDATER OSX 64 $revision $OSXDIR.tar.gz $OSXDIR.sh FDS/FDS6 test
