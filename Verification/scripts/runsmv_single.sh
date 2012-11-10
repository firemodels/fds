#!/bin/bash

function usage {
echo "runsmv_single.sh [-d -h -r -s size ]"
echo "Generates figures for Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of smokeview"
echo "-h - display this message"
echo "-r - run-time directory"
echo "-t - use test version of smokeview"
echo "-s size - use 32 or 64 bit (default) version of smokeview"
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
RUNOPTS=

while getopts 'dhtr:s:' OPTION
do
case $OPTION  in
  d)
   DEBUG=_db
   ;;
  h)
   usage;
   ;;
  r)
   RUNOPTS="$OPTARG"
   ;;
  t)
   TEST=_test
  ;;
  s)
   SIZE="$OPTARG"
   if [ $SIZE -eq 64 ] ; then
     SIZE=_64
   else
     SIZE=_32
   fi
  ;;
esac
done
shift $(($OPTIND-1))

if [ "$RUNOPTS" != "" ]; then
  RUNOPTS="-bindir $RUNOPTS"
fi
VERSION=$PLATFORM$TEST$SIZE$DEBUG
VERSION2=$PLATFORM$SIZE$DEBUG

export SVNROOT=~/FDS-SMV/

export SMV=$SVNROOT/SMV/Build/intel_$VERSION2/smokeview_$VERSION
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle"

export STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
export STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

source $STARTX
$SMV $RUNOPTS -runscript $1
echo command used: $SMV $RUNOPTS -redirect -runscript $1
source $STOPX
