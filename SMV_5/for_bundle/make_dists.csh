#!/bin/csh -f
set version=5.2.7_2602

set MACDIR=smv_$version\_osx
set LINUXDIR=smv_$version\_linux
set MACHOST=tiger.cfr.nist.gov

mkdir -p $MACDIR
cp readme.html $MACDIR/.
scp $MACHOST\:FDS-SMV/SMV_5/bin/smv5_osx $MACDIR/.
scp $MACHOST\:FDS-SMV/SMV_5/bin/smokezip_osx $MACDIR/.
rm -f $MACDIR.tar $MACDIR.tar.gz
tar cvf $MACDIR.tar $MACDIR/.
gzip $MACDIR.tar

mkdir -p $LINUXDIR
cp readme.html $LINUXDIR/.
cp ../bin/smv5_linux $LINUXDIR/.
cp ../bin/smokezip_linux $LINUXDIR/.
rm -f $LINUXDIR.tar $LINUXDIR.tar.gz
tar cvf $LINUXDIR.tar $LINUXDIR/.
gzip $LINUXDIR.tar
