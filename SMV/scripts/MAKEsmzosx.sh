#!/bin/bash
SVNROOT=~/$1

cd $SVNROOT/SMV/Build/smokezip/intel_osx_64
make -f ../Makefile clean >& /dev/null
./make_zip.sh
