#!/bin/bash
BUILD()
{
dir=$1
platform=intel_linux_64
script=$2
echo
echo "********** building $dir
echo
cd $dir/$platform
./$script
cd $BUILDDIR
}

CURDIR=`pwd`


cd ../../../../smv/Source
git clean -dxf

cd ../Build
git clean -dxf

BUILDDIR=`pwd`

BUILD LIBS makelibs.sh
BUILD background make_background.sh
BUILD dem2fds make_dem2fds.sh
BUILD fds2ascii make_fds2ascii.sh
BUILD hashfile make_hashfile.sh
BUILD smokediff make_smokediff.sh
BUILD wind2fds make_wind2fds.sh
BUILD smokeview make_smv.sh
