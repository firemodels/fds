#!/bin/bash

#PLATFORM=32
#IPLATFORM=ia32
PLATFORM=64
IPLATFORM=intel64
export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/INTEL_LINUX_$PLATFORM/smokeview_linux_$PLATFORM
export SMOKEZIP=$SVNROOT/Utilities/smokezip/INTEL_LINUX_$PLATFORM/smokezip_linux_$PLATFORM
export SMOKEDIFF=$SVNROOT/Utilities/smokediff/INTEL_LINUX_$PLATFORM/smokediff_linux_$PLATFORM
export BACKGROUND=$SVNROOT/Utilities/background/INTEL_LINUX_32/background
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
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
./SMV_Cases.sh
