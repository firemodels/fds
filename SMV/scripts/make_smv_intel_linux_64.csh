#!/bin/csh -f
set SVNROOT=~/FDS-SMV

#cd $SVNROOT/SMV/MACtiger2/sv5p0
cd $SVNROOT/SMV/Build/INTEL_LINUX_64
make -f  ../Makefile clean
./make_smv.csh
cd $SVNROOT/SMV/bin
