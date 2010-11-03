#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_OSX_32
make -f ../Makefile clean >& /dev/null
./make_zip.csh
cd $SVNROOT/Utilities/smokezip/INTEL_OSX_64
make -f ../Makefile clean >& /dev/null
./make_zip.csh
