#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/SMV/Build/wind2fds/intel_osx_64
make -f ../Makefile clean >& /dev/null
./make_wind.sh
