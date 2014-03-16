#!/bin/bash

function usage {
echo "runsmv_single.sh [-d -h -r -s size ]"
echo "Generates figures for Smokeview verification suite"
echo ""
echo "Options"
echo "-b - output build time in title"
echo "-d - use debug version of smokeview"
echo "-h - display this message"
echo "-m script - run smokeview script named script"
echo "-r - run-time directory"
echo "-t - use test version of smokeview"
echo "-s size - use 32 or 64 bit (default) version of smokeview"
echo "-x - invoke -volrender option in smokeview"
echo "-y startframe - invoke -startframe option in smokeview"
echo "-z skipframe - invoke -skipframe option in smokeview"
exit
}

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  PLATFORM=osx
else
  PLATFORM=linux
fi

SIZE=_64
DEBUG=
TEST=
SCRIPT=
RUNDIR=
RUNDIRECTORY=
SKIPFRAME=1
STARTFRAME=0
VOLRENDER=
TIME=

while getopts 'bdhm:r:s:txy:z:' OPTION
do
case $OPTION  in
  d)
   DEBUG=_db
   ;;
  h)
   usage;
   ;;
  m)
   SCRIPT="$OPTARG"
   ;;
  r)
   RUNDIR="$OPTARG"
   ;;
  b)
   TIME="-time"
   ;;
  t)
   TEST=_test
   ;;
  s)
   SIZE="$OPTARG"
   if [ $SIZE -eq 64 ]; then
     SIZE=_64
   else
     SIZE=_32
   fi
   ;;
  x)
   VOLRENDER=y
   ;;
  y)
   STARTFRAME="$OPTARG"
   ;;
  z)
   SKIPFRAME="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

if [ $# == 0 ] ; then
  usage;
  exit
fi

CASE=$1

# setup script

if [ "$SCRIPT" != "" ]; then
  SCRIPT="-script $SCRIPT"
  SCRIPTFILE=$SCRIPT
else
  SCRIPT="-script ${CASE}.ssf"
  SCRIPTFILE=${CASE}.ssf
fi

if [ "$RUNDIR" != "" ]; then
  RUNDIRECTORY=$RUNDIR
  RUNDIR="-bindir $RUNDIR"
fi

if [ "$VOLRENDER" == "y" ]; then
  VOLRENDER="-volrender"
  if [ "$STARTFRAME" != "" ]; then
    STARTFRAME="-startframe $STARTFRAME"
  fi
  if [ "$SKIPFRAME" != "" ]; then
    SKIPFRAME="-skipframe $SKIPFRAME"
  fi
else
  STARTFRAME=
  SKIPFRAME=
fi

VERSION=$PLATFORM$TEST$SIZE$DEBUG
VERSION2=$PLATFORM$SIZE$DEBUG

export SVNROOT=~/FDS-SMV/

export SMV=$SVNROOT/SMV/Build/intel_$VERSION2/smokeview_$VERSION
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle"
if [ "$RUNDIR" != "" ]; then
  SMVDIR="-bindir $RUNDIR"
fi

STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

echo setting up graphics environment
source $STARTX
echo "smokeview file: ${CASE}.smv"
echo "   script file: $SCRIPTFILE"
if [ "$RUNDIR" != "" ]; then
echo " bin directory: $RUNDIRECTORY"
fi
if [ "$VOLRENDER" == "y" ]; then
echo " *** volume rendering"
echo "   start frame: $STARTFRAME"
echo "    skip frame: $SKIPFRAME"
fi
$SMV $RUNDIR $SMVBINDIR $TIME $SCRIPT $VOLRENDER $STARTFRAME $SKIPFRAME $CASE >& /dev/null
echo shutting down graphics environment
source $STOPX
