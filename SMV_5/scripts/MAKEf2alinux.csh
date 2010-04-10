#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/fds2ascii/intel_linux_32
make_fds2ascii_linux32.csh
cd $SVNROOT/Utilities/fds2ascii/intel_linux_64
make_fds2ascii_linux64.csh
