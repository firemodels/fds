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
echo "-p path - specify path of the smokeview executable"
echo "-r - use release version of smokeview"
echo "-s size - use 32 or 64 bit (default) version of smokeview"
echo "-t - use test version of smokeview"
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
SMV_PATH=""

while getopts 'dhp:rs:t' OPTION
do
case $OPTION  in
  d)
   DEBUG=_dbg
   ;;
  h)
   usage;
   ;;
  p)
   SMV_PATH="$OPTARG"
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
  t)
   TEST=_test
  ;;
esac
done
shift $(($OPTIND-1))

export SVNROOT=`pwd`/..
if [ "$SMV_PATH" == "" ]; then
  SMV_PATH=$SVNROOT/SMV/Build/intel_$PLATFORM$SIZE
fi
export SMV=$SMV_PATH/smokeview_$PLATFORM$TEST$SIZE
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle/"
export BASEDIR=`pwd`

echo "erasing SCRIPT_FIGURES png files"
rm -f $SVNROOT/Manuals/FDS_Configuration_Management_Plan/SCRIPT_FIGURES/*.png
rm -f $SVNROOT/Manuals/FDS_Technical_Reference_Guide/SCRIPT_FIGURES/*.png
rm -f $SVNROOT/Manuals/FDS_User_Guide/SCRIPT_FIGURES/*.png
rm -f $SVNROOT/Manuals/FDS_Validation_Guide/SCRIPT_FIGURES/*.png
rm -f $SVNROOT/Manuals/FDS_Verificaiton_Guide/SCRIPT_FIGURES/*.png

source $SVNROOT/Utilities/Scripts/startXserver.sh
./FDS_Pictures.sh
source $SVNROOT/Utilities/Scripts/stopXserver.sh

echo FDS pictures created.

