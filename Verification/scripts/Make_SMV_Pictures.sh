#!/bin/bash

#PLATFORM=32
#IPLATFORM=ia32
PLATFORM=64
IPLATFORM=intel64
CURDIR=`pwd`
cd ..
export SVNROOT=`pwd`/..
#export SMV=$SVNROOT/SMV/Build/intel_linux__$PLATFORM/smokeview_linux_$PLATFORM
export SMV=$SVNROOT/SMV/Build/intel_linux_test_$PLATFORM/smokeview_linux_test_$PLATFORM
#export SMV=~/FDS/FDS5/bin/smv5_linux_64
export SMOKEZIP=$SVNROOT/Utilities/smokezip/intel_linux_$PLATFORM/smokezip_linux_$PLATFORM
export SMOKEDIFF=$SVNROOT/Utilities/smokediff/intel_linux_$PLATFORM/smokediff_linux_$PLATFORM
export BACKGROUND=$SVNROOT/Utilities/background/intel_linux_32/background
export RUNFDS=$SVNROOT/Utilities/Scripts/runsmv.sh
export BASEDIR=`pwd`

export SMVUG=$SVNROOT/Manuals/SMV_User_Guide
export SMVVG=$SVNROOT/Manuals/SMV_Verification_Guide

if ! [ -e $SMV ]; then
  echo "The file $SMV does not exist. Run aborted."
  exit
fi
if ! [ -e $SMOKEZIP ]; then
  echo "The file $SMOKEZIP does not exist. Run aborted."
  exit
fi
if ! [ -e $SMOKEDIFF ]; then
  echo "The file $SMOKEDIFF does not exist. Run aborted."
  exit
fi
if ! [ -e $BACKGROUND ]; then
  echo "The file $BACKGROUND does not exist. Run aborted."
  exit
fi


cd $SMVUG/SCRIPT_FIGURES
rm -f *.png
rm -f *.help
rm -f *.version
source ~/.bashrc_fds $IPLATFORM
$SMV -help > smokeview.help
$SMV -version > smokeview.version
$SMOKEZIP -help > smokezip.help
$SMOKEDIFF -help > smokediff.help
$SMOKEDIFF -v > smokediff.version
$BACKGROUND -help > background.help
$BACKGROUND -version > background.version

cd $SMVVG/SCRIPT_FIGURES
rm -f *.version
rm -f *.png
$SMV -version > smokeview.version

cd $SVNROOT/Verification
scripts/SMV_Cases.sh
cd $CURDIDR
