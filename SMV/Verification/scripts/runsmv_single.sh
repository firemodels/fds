#!/bin/bash

function usage {
echo "runsmv_single.sh [options ]"
echo "Run smokeview in batch mode"
echo ""
echo "Options"
echo "-b - output build time in title"
echo "-d - use debug version of smokeview"
echo "-h - display this message"
echo "-m script - run smokeview using the script named script"
echo "-r - smokeview bin directory"
echo "-t - use test version of smokeview"
echo "-s size - 32 or 64"
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
SKIPFRAME=1
STARTFRAME=0
VOLRENDER=
TIME=
BINDIR=

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
   BINDIR="$OPTARG"
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
VERSION2=$PLATFORM$SIZE

SVNROOT=~/FDS-SMV

SMOKEVIEW=$SVNROOT/SMV/Build/intel_$VERSION2/smokeview_$VERSION

if [ "$BINDIR" == "" ]; then
  BINDIR="$SVNROOT/SMV/for_bundle"
fi
SMVBINDIR="-bindir $BINDIR"

STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

source $STARTX
echo "     smokeview: $SMOKEVIEW"
echo "          case: ${CASE}.smv"
echo "        script: $SCRIPTFILE"
echo " bin directory: $BINDIR"
if [ "$VOLRENDER" == "y" ]; then
echo " *** volume rendering"
echo "   start frame: $STARTFRAME"
echo "    skip frame: $SKIPFRAME"
fi
$SMOKEVIEW $SMVBINDIR $TIME $SCRIPT $VOLRENDER $STARTFRAME $SKIPFRAME $CASE >& /dev/null
source $STOPX
