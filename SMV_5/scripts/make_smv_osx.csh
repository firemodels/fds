#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT
cd $SVNROOT/SMV_5/MACtiger2/sv5p0
make clean >& /dev/null
date >& $SVNROOT/SMV_5/bin/make_osx.out
make >>& $SVNROOT/SMV_5/bin/make_osx.out
cd $SVNROOT/SMV_5/bin
