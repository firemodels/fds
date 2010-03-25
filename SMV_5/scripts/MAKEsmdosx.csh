#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/INTEL_OSX_32
make -f ../Makefile clean >& /dev/null
./make_diff.csh
cd $SVNROOT/Utilities/smokediff/INTEL_OSX_64
make -f ../Makefile clean >& /dev/null
./make_diff.csh
