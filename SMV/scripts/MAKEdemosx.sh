#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/dem2fds/intel_osx_64
make -f ../Makefile clean >& /dev/null
./make_dem2fds.sh
