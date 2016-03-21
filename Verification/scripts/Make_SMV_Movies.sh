#!/bin/bash

size=_64
IB=
if [ "$FDSNETWORK" == "infiniband" ] ; then
IB=ib
fi

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx$size
else
  PLATFORM=linux$size
fi

CURDIR=`pwd`
cd ..
GITROOT=`pwd`
cd $CURDIR

export SMV=$GITROOT/SMV/Build/smokeview/intel_$PLATFORM/smokeview_$PLATFORM
FDSEXE=$GITROOT/FDS_Compilation/mpi_intel_$PLATFORM$IB/fds_mpi_intel_$PLATFORM$IB
RUNSMV="$GITROOT/Utilities/Scripts/runsmv.sh"
export SMVBINDIR="-bindir $GITROOT/SMV/for_bundle"
MAKEMOVIE=$GITROOT/Utilities/Scripts/make_movie.sh
STARTX=$GITROOT/Utilities/Scripts/startXserver.sh
STOPX=$GITROOT/Utilities/Scripts/stopXserver.sh
QFDS=$GITROOT/Utilities/Scripts/qfds.sh

export BASEDIR=`pwd`

if ! [ -e $SMV ]; then
  echo "*** Error: The program $SMV does not exist. Run aborted."
  exit
fi

underscore=_
mov=.m1v
VDIR=$GITROOT/Verification
INDIR=$GITROOT/Verification/Visualization/frames
WUIINDIR=$GITROOT/Verification/WUI/frames
OUTDIR=$GITROOT/Manuals/SMV_Summary/movies

rm -f $INDIR/*.png
rm -f $WUIINDIR/*.png
rm -f $OUTDIR/*.m1v
rm -f $OUTDIR/*.png

# make a movie
MKMOVIE()
{
  CASEDIR=$1
  BASEFILE=$2
  FRAMEDIR=$3

  cd $VDIR

# generate movie frames
  $RUNSMV -m -d $CASEDIR $BASEFILE

  cd $FRAMEDIR

# make movies out of frames generated above
  echo making $BASEFILE movie
  $MAKEMOVIE -o $OUTDIR ${BASEFILE}_movie  > /dev/null
}

# start background X server for picture generation
echo Starting background X server
source $STARTX

# create version string

cd $VDIR
$QFDS -e $FDSEXE -d Visualization -q terminal version2.fds

cd $VDIR
$RUNSMV -d Visualization version2

# The -m option assumes that a script
# named casename_movies.ssf exists for each 

# -------- plume5c movie -------------------

cd $VDIR

# generate movie frames
$RUNSMV -m -d Visualization plume5c

cd $INDIR

# make movies out of frames generated above
echo making plume5c_tslice movie
$MAKEMOVIE -o $OUTDIR -m plume5c_tslice plume5c_tslice > /dev/null

echo making plume5c_3dsmoke movie
$MAKEMOVIE -o $OUTDIR  -m plume5c_3dsmoke plume5c_3dsmoke  > /dev/null

echo making plume5c_vtslice movie
$MAKEMOVIE -o $OUTDIR -m plume5c_vtslice plume5c_vtslice > /dev/null

echo making plume5c_iso movie
$MAKEMOVIE -o $OUTDIR  -m plume5c_iso plume5c_iso  > /dev/null

echo making plume5c_tbound movie
$MAKEMOVIE -o $OUTDIR  -m plume5c_tbound plume5c_tbound > /dev/null

echo making plume5c_part movie
$MAKEMOVIE -o $OUTDIR  -m plume5c_part plume5c_part  > /dev/null

# -------- thouse5 movies -------------------

cd $VDIR

# generate movie frames
$RUNSMV -m -d Visualization thouse5

cd $INDIR

# make movies out of frames generated above
echo making thouse5_tslice movie
$MAKEMOVIE -o $OUTDIR  -m thouse5_tslice thouse5_tslice  > /dev/null

echo making thouse5_smoke3d movie
$MAKEMOVIE -o $OUTDIR  -m thouse5_smoke3d thouse5_smoke3d > /dev/null

# -------- BT10m_2x2km_LS movie -------------------

MKMOVIE WUI BT10m_2x2km_LS $WUIINDIR

# -------- hill_structure movie -------------------

MKMOVIE WUI hill_structure $WUIINDIR

# -------- levelset1 movie -------------------

MKMOVIE WUI levelset1 $WUIINDIR

# -------- wind_test1 movie -------------------

MKMOVIE WUI wind_test1 $WUIINDIR

# -------- tree_test2 movie -------------------

MKMOVIE WUI tree_test2 $WUIINDIR

source $STOPX
