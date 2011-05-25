#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_32
make -f ../Makefile clean >& /dev/null
./make_zip.sh
cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_64
make -f ../Makefile clean >& /dev/null
./make_zip.sh
