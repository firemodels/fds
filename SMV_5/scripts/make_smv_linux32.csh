#!/bin/csh -f
set SVNROOT=~/FDS-SMV
set rev=$1

cd $SVNROOT/SMV_5/source/smokeview
svn -r $rev update
cd $SVNROOT/SMV_5/Build/INTEL_LINUX_32
make -f ../Makefile clean >& /dev/null
date >& $SVNROOT/SMV_5/bin/make_intel_linux_32.out
make -f ../Makefile intel_linux_32 >>& $SVNROOT/SMV_5/bin/make_intel_linux_32.out
