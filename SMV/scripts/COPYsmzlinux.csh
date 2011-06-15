#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokezip/intel_linux_32
cp -v smokezip_linux_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokezip/intel_linux_64
cp -v smokezip_linux_64 $SVNROOT/SMV/bin/.
