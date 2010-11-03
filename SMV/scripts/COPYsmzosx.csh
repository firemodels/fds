#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_OSX_32
cp -v smokezip_osx_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokezip/INTEL_OSX_64
cp -v smokezip_osx_64 $SVNROOT/SMV/bin/.
