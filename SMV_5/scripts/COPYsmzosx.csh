#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_OSX_32
cp smokezip_osx_32 $SVNROOT/SMV_5/bin/.
cd $SVNROOT/Utilities/smokezip/INTEL_OSX_64
cp smokezip_osx_64 $SVNROOT/SMV_5/bin/.
