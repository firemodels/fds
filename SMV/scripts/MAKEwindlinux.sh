#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/wind2fds/intel_linux_64
make -f ../Makefile clean >& /dev/null
./make_wind.sh
