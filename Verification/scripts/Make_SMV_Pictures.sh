#!/bin/bash

function usage {
echo "Make_SMV_Pictures.sh [-d -h -r -s size ]"
echo "Generates figures for Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of smokeview"
echo "-h - display this message"
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

while getopts 'dhts:' OPTION
do
case $OPTION  in
  d)
   DEBUG=_db
   ;;
  h)
   usage;
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


VERSION=$PLATFORM$TEST$SIZE$DEBUG
VERSION2=$PLATFORM$SIZE
IPLATFORM=intel64
CURDIR=`pwd`
cd ../..
export SVNROOT=`pwd`
cd $CURDIR/..

export SMV=$SVNROOT/SMV/Build/intel_$VERSION2$DEBUG/smokeview_$VERSION
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle"

export SMOKEZIP=$SVNROOT/Utilities/smokezip/intel_$VERSION2/smokezip_$VERSION2
export SMOKEDIFF=$SVNROOT/Utilities/smokediff/intel_$VERSION2/smokediff_$VERSION2
export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM\_32/background
export STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
export STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

echo Generating smokeview images using:
echo smokeview : $SMV $SMVBINDIR
echo smokezip  : $SMOKEZIP
echo smokediff : $SMOKEDIFF
echo background: $BACKGROUND
echo

export RUNFDS=$SVNROOT/Utilities/Scripts/runsmv.sh
export RUNCFAST=$SVNROOT/Utilities/Scripts/runsmv.sh
export BASEDIR=`pwd`

export SMVUG=$SVNROOT/Manuals/SMV_User_Guide
export SMVVG=$SVNROOT/Manuals/SMV_Verification_Guide
SUMMARY=$SVNROOT/Manuals/Verification_Summary

if ! [ -e $SMV ]; then
  echo "The file $SMV does not exist. Run aborted."
  exit
fi
if ! [ -e $SMOKEZIP ]; then
  echo "The file $SMOKEZIP does not exist. Run aborted."
  exit
fi
if ! [ -e $SMOKEDIFF ]; then
  echo "The file $SMOKEDIFF does not exist. Run aborted."
  exit
fi
if ! [ -e $BACKGROUND ]; then
  echo "The file $BACKGROUND does not exist. Run aborted."
  exit
fi

cd $SMVUG/SCRIPT_FIGURES
rm -f *.png
rm -f *.help
rm -f smokeview.version
rm -f smokediff.version
rm -f smokezip.version
rm -f background.version
rm -f $SUMMARY/images/*.png
source ~/.bashrc_fds $IPLATFORM

$SMV -help > smokeview.help
$SMOKEZIP -help > smokezip.help
$SMOKEDIFF -help > smokediff.help
$BACKGROUND -help > background.help

$SMV -version > smokeview.version
$SMOKEZIP -v > smokezip.version
$SMOKEDIFF -v > smokediff.version
$BACKGROUND -version > background.version

cd $SMVVG/SCRIPT_FIGURES
rm -f *.version
rm -f *.png
$SMV -version > smokeview.version

cd $SVNROOT/Verification/Visualization
echo Converting particles to isosurfaces in case plumeiso
$SMOKEZIP -r -part2iso plumeiso

cd $SVNROOT/Verification/Visualization
echo Differencing cases plume5c and plume5cdelta
$SMOKEDIFF -r plume5c plume5cdelta
echo Differencing cases thouse5 and thouse5delta
$SMOKEDIFF -r thouse5 thouse5delta

echo Generating images

source $STARTX
cd $SVNROOT/Verification
scripts/SMV_Cases.sh

cd $SVNROOT/Verification
scripts/SMV_DIFF_Cases.sh
cd $CURDIDR
source $STOPX

# copy generated images to web summary directory

cp $SMVUG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
cp $SMVVG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
