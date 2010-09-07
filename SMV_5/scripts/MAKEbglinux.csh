#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/background/INTEL_LINUX_32
make -f ../Makefile clean >& /dev/null
./make_background.csh
