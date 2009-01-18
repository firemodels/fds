#!/bin/csh -f
set SVNROOT=~/FDS-SMV

#cd $SVNROOT/SMV_5/MACtiger2/sv5p0
cd $SVNROOT/SMV_5/Build/INTEL_OSX_32
make -f ../Makefile clean >& /dev/null
date >& $SVNROOT/SMV_5/bin/make_intel_osx_32.out
./make_smv.csh >>& $SVNROOT/SMV_5/bin/make_intel_osx_32.out
cd $SVNROOT/SMV_5/bin
