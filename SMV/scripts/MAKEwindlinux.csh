#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/wind2fds/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_wind.sh
