#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/intel_osx_64
make -f ../Makefile clean >& /dev/null
./make_zip.sh
