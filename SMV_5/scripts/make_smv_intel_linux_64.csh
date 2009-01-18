#!/bin/csh -f
set SVNROOT=~/FDS-SMV

#cd $SVNROOT/SMV_5/MACtiger2/sv5p0
cd $SVNROOT/SMV_5/Build/INTEL_LINUX_64
make -f  ../Makefile clean >& /dev/null
date >& $SVNROOT/SMV_5/bin/make_intel_linux_64.out
./make_smv.csh >>& $SVNROOT/SMV_5/bin/make_intel_linux_64.out
cd $SVNROOT/SMV_5/bin
