#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/INTEL_OSX_32
cp -v smokediff_osx_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokediff/INTEL_OSX_64
cp -v smokediff_osx_64 $SVNROOT/SMV/bin/.
