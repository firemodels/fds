#!/bin/bash
SVNROOT=~/$1
echo arg1=$1
cd $SVNROOT/SMV/Build/background/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_background.sh
