#!/bin/bash
curdir=`pwd`
testimage=plume5c_test.png
rm -f $testimage
cd ../..
reporoot=`pwd`
visdir=$reporoot/Verification/Visualization
if [ "`uname`" == "Darwin" ]; then
  SMV="$reporoot/Utilities/Scripts/smokeview.sh -e $reporoot/SMV/Build/intel_osx_64/smokeview_osx_64"
else
  SMV="$reporoot/Utilities/Scripts/smokeview.sh -e $reporoot/SMV/Build/intel_linux_64/smokeview_linux_64"
fi
cd $visdir
$SMV -s plume5c2.ssf plume5c
cd $curdir
if [ -e $testimage ]; then
  echo test image creation succeeded
else
  echo test image creation failed
fi
cd $curdir
