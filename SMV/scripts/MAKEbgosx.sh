#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/background/intel_osx_64
make -f ../Makefile clean >& /dev/null
./make_background.sh
