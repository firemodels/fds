#!/bin/csh -f
set SVNROOT=~/FDS-SMV

cd $SVNROOT/Utilities/smokediff/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_diff.sh
