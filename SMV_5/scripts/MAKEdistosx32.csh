#!/bin/csh -f
set version=$1
set revision=$2
set SVNROOT=~/FDS-SMV
set OSXHOST=bluesky.cfr.nist.gov

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set OSXDIR=smv_$version\_osx32

cd $FORBUNDLE

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp readme.html $OSXDIR/Documentation/release_notes.html

cp $FORBUNDLE/objects.svo $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokediff_osx_32 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smv5_osx_32 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokezip_osx_32 $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar
