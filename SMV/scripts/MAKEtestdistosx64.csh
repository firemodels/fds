#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV
set REMOTESVNROOT=FDS-SMV
set OSXHOST=$2

set SMVDIR=$REMOTESVNROOT/SMV/Build/INTEL_OSX_TEST_64
set SMZDIR=$REMOTESVNROOT/Utilities/smokezip/INTEL_OSX_64
set SMDDIR=$REMOTESVNROOT/Utilities/smokediff/INTEL_OSX_64
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set OSXDIR=smv_test\_$revision\_osx_64

cd $FORBUNDLE

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp note.txt $OSXDIR/Documentation/.

cp $FORBUNDLE/objects.svo $OSXDIR/.
scp $OSXHOST\:$SMVDIR/smokeview_osx_test_64 $OSXDIR/.
scp $OSXHOST\:$SMZDIR/smokezip_osx_64 $OSXDIR/.
scp $OSXHOST\:$SMDDIR/smokediff_osx_64 $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar
