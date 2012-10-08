#!/bin/bash

PLATFORM=64
export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/intel_linux_$PLATFORM/smokeview_linux_$PLATFORM
export RUNSMV="$SVNROOT/Utilities/Scripts/runsmv.sh"
export MAKEMOVIE=$SVNROOT/Utilities/Scripts/makemovie.sh
export STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
export STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

export BASEDIR=`pwd`

if ! [ -e $SMV ]; then
  echo "*** Error: The program $SMV does not exist. Run aborted."
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

# start background X server for picture generation
echo Starting background X server
source $STARTX

# This script assume that a smokeview scriptfile 
# named casename_movies.ssf exists for each 
$RUNSMV -m Visualization plume5c

cd $INDIR

echo making plume5c_tslice movie
$MAKEMOVIE plume5c_tslice $OUTDIR > /dev/null
echo making plume5c_3dsmoke movie
$MAKEMOVIE plume5c_3dsmoke $OUTDIR > /dev/null
echo making plume5c_vtslice movie
$MAKEMOVIE plume5c_vtslice $OUTDIR > /dev/null
echo making plume5c_iso movie
$MAKEMOVIE plume5c_iso $OUTDIR > /dev/null
echo making plume5c_tbound movie
$MAKEMOVIE plume5c_tbound $OUTDIR > /dev/null
echo making plume5c_part movie
$MAKEMOVIE plume5c_part $OUTDIR > /dev/null

cd $VDIR

$RUNSMV -m Visualization thouse5

cd $INDIR

echo making thouse5_tslice movie
$MAKEMOVIE thouse5_tslice $OUTDIR > /dev/null
echo making thouse5_smoke3d movie
$MAKEMOVIE thouse5_smoke3d $OUTDIR > /dev/null

source $STOPX
