#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT
cd $SVNROOT/SMV_5/MACtiger2/sv5p0
make clean >& /dev/null
make >& $SVNROOT/SMV_5/bin/make_osx.out
cd $SVNROOT/SMV_5/bin
ls -l smv5_osx
