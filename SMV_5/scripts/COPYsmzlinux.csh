#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_32
cp smokezip_linux_32 $SVNROOT/SMV_5/bin/.
cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_64
cp smokezip_linux_64 $SVNROOT/SMV_5/bin/.
