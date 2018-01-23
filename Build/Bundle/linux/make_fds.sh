#!/bin/bash
CURDIR=`pwd`
cd ../../../Build/impi_intel_linux_64
git clean -dxf
echo
echo "********** build FDS"
echo
./make_fds.sh
cd $CURDIR
