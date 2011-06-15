#!/bin/csh -f
set version=$1
set revision=$2
set SVNROOT=FDS-SMV
set OSXHOST=$3

set SMOKEVIEWDIR=$SVNROOT/SMV/bin
set SMOKEZIPDIR=$SVNROOT/Utilities/smokezip/intel_osx_32
set SMOKEDIFFDIR=$SVNROOT/Utilities/smokediff/intel_osx_32
set FORBUNDLE=~/$SVNROOT/SMV/for_bundle
set OSXDIR=smv_$version\_osx32

cd $FORBUNDLE

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp readme.html $OSXDIR/Documentation/release_notes.html

cp $FORBUNDLE/objects.svo $OSXDIR/.
scp $OSXHOST\:$SMOKEDIFFDIR/smokediff_osx_32 $OSXDIR/.
scp $OSXHOST\:$SMOKEVIEWDIR/smv5_osx_32 $OSXDIR/.
scp $OSXHOST\:$SMOKEZIPDIR/smokezip_osx_32 $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar
