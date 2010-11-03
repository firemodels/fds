#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_32
cp -v smokezip_linux_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokezip/INTEL_LINUX_64
cp -v smokezip_linux_64 $SVNROOT/SMV/bin/.
