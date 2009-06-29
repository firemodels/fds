#!/bin/csh -f
set version=$1
set SVNROOT=~/FDS-SMV
set OSXHOST=tiger.cfr.nist.gov

set BINDIR=$SVNROOT/SMV_5/bin
set FORBUNDLE=$SVNROOT/SMV_5/for_bundle
set OSXDIR=smv_$version\_osx
set LINUXDIR=smv_$version\_linux
set LINUXDIR64=smv_$version\_linux_64

cd $FORBUNDLE

mkdir -p $OSXDIR
mkdir -p $OSXDIR/Documentation
cp readme.html $OSXDIR/Documentation/.

scp $OSXHOST\:FDS-SMV/SMV_5/bin/smv5_osx_32 $OSXDIR/.
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokezip_osx $OSXDIR/.
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar

mkdir -p $LINUXDIR
mkdir -p $LINUXDIR/Documentation
cp readme.html $LINUXDIR/Documentation/.
cp $BINDIR/smv5_linux_32 $LINUXDIR/.
cp $BINDIR/smokezip_linux $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar

rm -rf $LINUXDIR64
mkdir -p $LINUXDIR64
mkdir -p $LINUXDIR64/Documentation
cp readme.html $LINUXDIR64/Documentation/.

cp $BINDIR/smv5_linux_64 $LINUXDIR64/.
cp $BINDIR/smokezip_linux_64 $LINUXDIR64/.
rm -f $LINUXDIR64.tar $LINUXDIR64.tar.gz
tar cvf $LINUXDIR64.tar $LINUXDIR64/.
gzip $LINUXDIR64.tar
