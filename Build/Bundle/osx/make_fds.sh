#!/bin/bash
CURDIR=`pwd`
cd ../../../Build/mpi_intel_osx_64
git clean -dxf
echo
echo "********** build FDS"
echo
./make_fds.sh
cd $CURDIR
