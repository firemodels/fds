#!/bin/csh -f
set revision=$1
set SVNROOT=~/FDS-SMV
set OSXHOST=bluesky.cfr.nist.gov

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set OSXDIR=smv_test\_$revision\_osx_64

cd $FORBUNDLE

rm -rf $OSXDIR
mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp note.txt $OSXDIR/Documentation/.

cp $FORBUNDLE/objects.svo $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smv5_osx_test_64 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokezip_osx_64 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokediff_osx_64 $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar
