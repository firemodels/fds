#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/intel_linux_32
cp -v smokediff_linux_32 $SVNROOT/SMV/bin/.
cd $SVNROOT/Utilities/smokediff/intel_linux_32
cp -v smokediff_linux_64 $SVNROOT/SMV/bin/.
