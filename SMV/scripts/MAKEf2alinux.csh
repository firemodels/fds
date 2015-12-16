#!/bin/csh -f
set SVNROOT=~/$1

cd $SVNROOT/Utilities/fds2ascii/intel_linux_64
./make_fds2ascii.sh
