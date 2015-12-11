#!/bin/bash
SVNROOT=~/$1
echo arg1=$1
echo $SVNROOT=$SVNROOT
cd $SVNROOT/Utilities/background/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_background.sh
