#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/fds2ascii/intel_osx_32
./make_fds2ascii_osx32.csh
cd $SVNROOT/Utilities/fds2ascii/intel_osx_64
./make_fds2ascii_osx64.csh
