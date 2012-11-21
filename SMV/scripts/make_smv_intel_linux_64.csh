#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/SMV/Build/intel_linux_64
make -f  ../Makefile clean
./make_smv.csh
cd $SVNROOT/SMV/bin
