#!/bin/bash -f

# This script generates SMV pictures from the
# FDS Verification Cases on a Linux or OS X machine

function usage {
echo "Make_SMV_Pictures.sh [-d -h -r -s size ]"
echo "Generates Smokeview figures from FDS verification suite"
echo ""
echo "Options"
echo "-d - use debug version of smokeview"
echo "-h - display this message"
echo "-r - use release version of smokeview"
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
TEST=_test

while getopts 'dhrs:' OPTION
do
case $OPTION  in
  d)
   DEBUG=_dbg
   ;;
  h)
   usage;
   ;;
  r)
   TEST=
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

VERSION=$PLATFORM$TEST$SIZE$DEBUG

export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/intel_$VERSION/smokeview_$VERSION
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle/"
export BASEDIR=`pwd`

./FDS_Pictures.sh

echo FDS cases submitted

