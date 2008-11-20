#!/bin/csh -f
set version=5.3_2700
set SVNROOT=~/FDS-SMV
set MACHOST=tiger.cfr.nist.gov

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set MACDIR=smv_$version\_osx
set LINUXDIR=smv_$version\_linux

cd $FORBUNDLE

mkdir -p $MACDIR
cp readme.html $MACDIR/.
scp $MACHOST\:FDS-SMV/SMV_5/bin/smv5_osx $MACDIR/.
scp $MACHOST\:FDS-SMV/SMV_5/bin/smokezip_osx $MACDIR/.
rm -f $MACDIR.tar $MACDIR.tar.gz
tar cvf $MACDIR.tar $MACDIR/.
gzip $MACDIR.tar

mkdir -p $LINUXDIR
cp readme.html $LINUXDIR/.
cp $BINDIR/smv5_linux $LINUXDIR/.
cp $BINDIR/smokezip_linux $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
