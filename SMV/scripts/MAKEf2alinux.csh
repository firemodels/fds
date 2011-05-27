#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/fds2ascii/intel_linux_32
./make_fds2ascii.sh
cd $SVNROOT/Utilities/fds2ascii/intel_linux_64
./make_fds2ascii.sh
