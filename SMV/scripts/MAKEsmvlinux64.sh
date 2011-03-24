#!/bin/bash
SVNROOT=~/FDS-SMV
rev=$1

cd $SVNROOT/SMV/source/smokeview
svn -r $rev update
cd $SVNROOT/SMV/Build/INTEL_LINUX_64
make -f ../Makefile clean >& /dev/null
./make_smv.sh
