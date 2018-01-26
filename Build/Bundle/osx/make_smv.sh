#!/bin/bash
BUILDFDS()
{
dir=$1
platform=intel_osx_64
script=make_${dir}.sh

echo
echo "********** building $dir"
echo
cd $CURDIR/../../../Utilities
cd $dir/$platform
./$script
cd $BUILDDIR
}

BUILD()
{
dir=$1
platform=intel_osx_64
script=make_${dir}.sh

echo
echo "********** building $dir"
echo
cd $dir/$platform
./$script
cd $BUILDDIR
}

CURDIR=`pwd`
bundlelog=bundle.log

cd ../../../../smv/Source
git clean -dxf >& /dev/null

cd ../Build
git clean -dxf >& /dev/null

BUILDDIR=`pwd`
bundlelog=$CURDIR/bundle.log
echo > $bundlelog

BUILD LIBS         | tee -a $bundlelog
BUILD background   | tee -a $bundlelog
BUILD dem2fds      | tee -a $bundlelog
BUILDFDS fds2ascii | tee -a $bundlelog
BUILD hashfile     | tee -a $bundlelog
BUILD smokediff    | tee -a $bundlelog
BUILD smokezip     | tee -a $bundlelog
BUILD wind2fds     | tee -a $bundlelog
BUILD smokeview    | tee -a $bundlelog
echo "smokeview builds complete" | tee -a $bundlelog
