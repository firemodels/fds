#!/bin/bash
testimage=plume5c_test.png
curdir=`pwd`
rm -f $testimage
cd ../..
reporoot=`pwd`
visdir=$reporoot/Verification/Visualization
SMV=$repo/Utilties/runsmv.sh
cd $visdir
$SMV -s plume5c2.ssf plume5c
cd $curdir
if [ -e $testimage ]; then
  echo test image creation succeeded
else
  echo test image creation failed
fi

