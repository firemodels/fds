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

underscore=_
mov=.m1v


cd $SVNROOT/Verification

# This script assume that a smokeview scriptfile 
# named casename_movies.ssf exists for each 

# plume5c movies

$RUNSMV Visualization plume5c
cd $SVNROOT/Verification/Visualization
movie=plume5c_tslice
$MAKEMOVIE $movie$underscore $movie
mv $movie$mov ../../Manuals/Movies/.

cd $SVNROOT/Verification
