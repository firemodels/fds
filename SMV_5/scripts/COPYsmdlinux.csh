#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/INTEL_LINUX_32
cp -v smokediff_linux_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokediff/INTEL_LINUX_64
cp -v smokediff_linux_64 $SVNROOT/SMV/bin/.
