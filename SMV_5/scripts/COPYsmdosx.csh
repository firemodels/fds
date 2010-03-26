#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/INTEL_OSX_32
cp smokediff_osx_32 $SVNROOT/SMV_5/bin/.
cd $SVNROOT/Utilities/smokediff/INTEL_OSX_64
cp smokediff_osx_64 $SVNROOT/SMV_5/bin/.
