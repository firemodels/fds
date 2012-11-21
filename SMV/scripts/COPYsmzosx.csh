#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/intel_osx_32
cp -v smokezip_osx_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokezip/intel_osx_64
cp -v smokezip_osx_64 $SVNROOT/SMV/bin/.
