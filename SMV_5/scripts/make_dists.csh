#!/bin/csh -f
set version=$1
set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set OSXDIR=smv_$version\_osx
set LINUXDIR=smv_$version\_linux

cd $FORBUNDLE

mkdir -p $OSXDIR
cp readme.html $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smv5_osx $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokezip_osx $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar

mkdir -p $LINUXDIR
cp readme.html $LINUXDIR/.
cp $BINDIR/smv5_linux $LINUXDIR/.
cp $BINDIR/smokezip_linux $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
