#!/bin/bash

#PLATFORM=32
PLATFORM=64
export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/intel_linux_$PLATFORM/smokeview_linux_$PLATFORM
#export SMV=~/FDS/FDS6/bin/smokeview_linux_$PLATFORM
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
export OUTDIR=$SVNROOT/Manuals/SMV_Animations

rm -f $INDIR/*.png
rm -f $OUTDIR/*.m1v
cd $VDIR

# This script assume that a smokeview scriptfile 
# named casename_movies.ssf exists for each 

# plume5c movies

$RUNSMV Visualization plume5c
cd $INDIR

$MAKEMOVIE plume5c_tslice $OUTDIR
$MAKEMOVIE plume5c_3dsmoke $OUTDIR
$MAKEMOVIE plume5c_vtslice $OUTDIR
$MAKEMOVIE plume5c_iso $OUTDIR
$MAKEMOVIE plume5c_tbound $OUTDIR
$MAKEMOVIE plume5c_part $OUTDIR

cd $VDIR
