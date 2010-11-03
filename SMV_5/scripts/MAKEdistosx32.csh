#!/bin/csh -f
set version=$1
set revision=$2
set SVNROOT=~/FDS-SMV
set OSXHOST=$3

set BINDIR=$SVNROOT/SMV/bin
set FORBUNDLE=$SVNROOT/SMV/for_bundle
set OSXDIR=smv_$version\_osx32

cd $FORBUNDLE

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp readme.html $OSXDIR/Documentation/release_notes.html

cp $FORBUNDLE/objects.svo $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV/bin/smokediff_osx_32 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV/bin/smv5_osx_32 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV/bin/smokezip_osx_32 $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar
