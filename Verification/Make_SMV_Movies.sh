#!/bin/bash

#PLATFORM=32
PLATFORM=64
export SVNROOT=`pwd`/..
#export SMV=$SVNROOT/SMV/Build/INTEL_LINUX_$PLATFORM/smokeview_linux_$PLATFORM
export SMV=~/FDS/FDS6/bin/smokeview_linux_$PLATFORM
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv_movie.sh
export MAKEMOVIE=$SVNROOT/Utilities/Scripts/makemovie.sh
export BASEDIR=`pwd`

if ! [ -e $SMV ]; then
  echo "The file $SMV does not exist. Run aborted."
  exit
fi

underscore=_
mov=.m1v
VDIR=$SVNROOT/Verification
INDIR=$SVNROOT/Verification/Visualization/frames
OUTDIR=$SVNROOT/Manuals/SMV_Animations

rm -f $INDIR/*.png
rm -f $OUTDIR/*.m1v
cd $VDIR

# This script assume that a smokeview scriptfile 
# named casename_movies.ssf exists for each 

# plume5c movies

$RUNSMV Visualization plume5c
cd $INDIR

movie=plume5c_tslice
$MAKEMOVIE $movie$underscore $movie
mv $movie$mov $OUTDIR/.

movie=plume5c_3dsmoke
$MAKEMOVIE $movie$underscore $movie
mv $movie$mov $OUTDIR/.


cd $VDIR
