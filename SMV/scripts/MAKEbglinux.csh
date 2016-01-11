#!/bin/csh -f
set SVNROOT=~/$1

cd $SVNROOT/Utilities/background/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_background.sh
