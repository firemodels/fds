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

rm -rf $OSXDIR
mkdir -p $OSXDIR

cp readme.html $OSXDIR/release_notes.html
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smv5_osx_32 $OSXDIR/smv5_osx
scp $OSXHOST\:FDS-SMV/SMV_5/bin/smokezip_osx $OSXDIR/.
scp $OSXHOST\:FDS-SMV/Utilities/smokediff/INTEL_OSX_32/smokediff_osx_32 $OSXDIR/smokediff_osx
rm -f $OSXDIR.tar $OSXDIR.tar.gz
tar cvf $OSXDIR.tar $OSXDIR/.
gzip $OSXDIR.tar

rm -rf $LINUXDIR
mkdir -p $LINUXDIR

cp readme.html $LINUXDIR/release_notes.html
cp $BINDIR/smv5_linux_32 $LINUXDIR/smv5_linux
cp $BINDIR/smokezip_linux $LINUXDIR/.
cp $SVNROOT/Utilities/smokediff/INTEL_LINUX_32/smokediff_linux_32 $LINUXDIR/smokediff_linux
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar

rm -rf $LINUXDIR64
mkdir -p $LINUXDIR64

cp readme.html $LINUXDIR64/release_notes.html
cp $BINDIR/smv5_linux_64 $LINUXDIR64/.
cp $BINDIR/smokezip_linux_64 $LINUXDIR64/.
cp $SVNROOT/Utilities/smokediff/INTEL_LINUX_64/smokediff_linux_64 $LINUXDIR64/smokediff_linux_64
rm -f $LINUXDIR64.tar $LINUXDIR64.tar.gz
tar cvf $LINUXDIR64.tar $LINUXDIR64/.
gzip $LINUXDIR64.tar
