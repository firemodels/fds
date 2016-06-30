#!/bin/bash
curdir=`pwd`
testimage=plume5c_test.png
rm -f $testimage
cd ../..
# define location of various directories
reporoot=`pwd`
visdir=$reporoot/Verification/Visualization
if [ "`uname`" == "Darwin" ]; then
  SMVDIR=$reporoot/SMV/Build/intel_osx_64
  SMVPROG=smokeview_osx_64
  SMV="$reporoot/SMV/Utilities/Scripts/smokeview.sh -e $SMVDIR/$SMVPROG"
else
  SMVDIR=$reporoot/SMV/Build/intel_linux_64
  SMVPROG=smokeview_linux_64
  SMV="$reporoot/SMV/Utilities/Scripts/smokeview.sh -e $SMVDIR/$SMVPROG"
fi

# building smokeview if it does not exist
if [ ! -e $SMVDIR/$SMVPROG ]; then
  echo smokeview does not exist - building smokeview
  cd $SMVDIR
  ./make_smv.sh
fi

cd $visdir
# make sure case exists
if [ ! -e plume5c.smv ]; then
  echo plume5c.smv does not exist
  echo run the case plume5c.fds 
  exit
fi

# generate image
$SMV -s plume5c2.ssf plume5c
cd $curdir

# check to see if image exists
if [ -e $testimage ]; then
  echo ""
  echo "*** Success! The image $testimage was created."
else
  echo ""
  echo "*** Failure.  The image $testimage was not created."
fi
cd $curdir
