#!/bin/csh -f
set SVNROOT=~/FDS-SMV
set rev=$1

cd $SVNROOT/SMV_5/source/Smokeview
svn -r $rev update
cd $SVNROOT/SMV_5/Build/INTEL_LINUX_64
make -f ../Makefile clean >& /dev/null
date >& $SVNROOT/SMV_5/bin/make_intel_linux_64.out
make -f ../Makefile intel_linux_64 >>& $SVNROOT/SMV_5/bin/make_intel_linux_64.out
