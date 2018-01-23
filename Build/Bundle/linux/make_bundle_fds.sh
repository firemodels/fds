#!/bin/bash
CURDIR=`pwd`
cd ../../../Build/impi_intel_linux_64
git clean -dxf
./make_fds.sh
cd $CURDIR
