#!/bin/bash

#PLATFORM=32
#IPLATFORM=ia32
PLATFORM=64
IPLATFORM=intel64
export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/INTEL_LINUX_$PLATFORM/smokeview_linux_$PLATFORM
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv_movie.sh
export MAKEMOVIE=$SVNROOT/Utilities/Scripts/makemovie.sh
export BASEDIR=`pwd`

if ! [ -e $SMV ]; then
  echo "The file $SMV does not exist. Run aborted."
  exit
fi


cd $SVNROOT/Verification

#$RUNSMV Visualization plume5c
cd $SVNROOT/Verification/Visualization
source ~/.bashrc_fds ia32
$MAKEMOVIE plume5c_tslice plume5c_tslice

cd $SVNROOT/Verification
