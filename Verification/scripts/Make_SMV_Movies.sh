#!/bin/bash

size=64

while getopts 'p:' OPTION
do
case $OPTION in
  p)
   size="$OPTARG"
   ;;
esac
#shift
done

if [ "$size" != "32" ]; then
  size=64
fi
size=_$size

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
else
  PLATFORM=linux$size
fi

export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/intel_$PLATFORM/smokeview_$PLATFORM
export RUNSMV="$SVNROOT/Utilities/Scripts/runsmv.sh"
export MAKEMOVIE=$SVNROOT/Utilities/Scripts/make_movie.sh
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

# The -m option assumes that a script
# named casename_movies.ssf exists for each 

# generate movie frames
$RUNSMV -m Visualization plume5c

cd $INDIR

# make movies out of frames generated above
echo making plume5c_tslice movie
$MAKEMOVIE -o $OUTDIR plume5c_tslice  > /dev/null
echo making plume5c_3dsmoke movie
$MAKEMOVIE -o $OUTDIR plume5c_3dsmoke  > /dev/null
echo making plume5c_vtslice movie
$MAKEMOVIE -o $OUTDIR plume5c_vtslice > /dev/null
echo making plume5c_iso movie
$MAKEMOVIE -o $OUTDIR plume5c_iso  > /dev/null
echo making plume5c_tbound movie
$MAKEMOVIE -o $OUTDIR plume5c_tbound  > /dev/null
echo making plume5c_part movie
$MAKEMOVIE -o $OUTDIR plume5c_part  > /dev/null

cd $VDIR

# generate movie frames
$RUNSMV -m Visualization thouse5

cd $INDIR

# make movies out of frames generated above
echo making thouse5_tslice movie
$MAKEMOVIE -o $OUTDIR thouse5_tslice  > /dev/null
echo making thouse5_smoke3d movie
$MAKEMOVIE -o $OUTDIR thouse5_smoke3d  > /dev/null

source $STOPX
