#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/SMV/Build/smokediff/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_diff.sh
