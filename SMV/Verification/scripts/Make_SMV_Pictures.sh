#!/bin/bash

function usage {
echo "Make_SMV_Pictures.sh [-d -h -r -s size ]"
echo "Generates figures for Smokeview verification suite"
echo ""
echo "Options"
echo "-d - use debug version of smokeview"
echo "-g - only generate geometry case images"
echo "-h - display this message"
echo "-i - use installed version of smokeview"
echo "-t - use test version of smokeview"
echo "-s size - use 32 or 64 bit (default) version of smokeview"
echo "-W - only generate WUI case images"
exit
}

is_file_installed()
{
  program=$1
  notfound=`$program -help 2>&1 | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    echo "The program $SMV is not available. Run aborted."
    exit
  fi
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
use_installed="0"
RUN_SMV=1
RUN_GEOM=0
RUN_WUI=1

while getopts 'dghis:tW' OPTION
do
case $OPTION  in
  d)
   DEBUG=_db
   ;;
  g)
   RUN_SMV=0
   RUN_GEOM=1
   RUN_WUI=0
   ;;
  h)
   usage;
   ;;
  i)
   use_installed="1"
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
  W)
   RUN_SMV=0
   RUN_GEOM=0
   RUN_WUI=1
   ;;
esac
done
shift $(($OPTIND-1))


VERSION=$PLATFORM$TEST$SIZE$DEBUG
VERSION2=$PLATFORM$SIZE
CURDIR=`pwd`
cd ../../..
export SVNROOT=`pwd`
cd $CURDIR/..

if [ "$use_installed" == "1" ] ; then
  export SMV=smokeview
  export SMOKEZIP=smokediff
  export SMOKEDIFF=smokediff
  export WIND2FDS=wind2fds
  export BACKGROUND=background
else
  export SMV=$SVNROOT/SMV/Build/smokeview/intel_$VERSION2/smokeview_$VERSION
  export SMOKEZIP=$SVNROOT/SMV/Build/smokezip/intel_$VERSION2/smokezip_$VERSION2
  export SMOKEDIFF=$SVNROOT/SMV/Build/smokediff/intel_$VERSION2/smokediff_$VERSION2
  export WIND2FDS=$SVNROOT/SMV/Build/wind2fds/intel_$VERSION2/wind2fds_$VERSION2
  export BACKGROUND=$SVNROOT/SMV/Build/background/intel_$VERSION2/background
fi

export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle"

export STARTX=$SVNROOT/Utilities/Scripts/startXserver.sh
export STOPX=$SVNROOT/Utilities/Scripts/stopXserver.sh

echo Generating smokeview images using:
echo smokeview : $SMV $SMVBINDIR
echo smokezip  : $SMOKEZIP
echo smokediff : $SMOKEDIFF
echo background: $BACKGROUND
echo

RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export QFDS=$RUNSMV
export RUNCFAST=$RUNSMV
export BASEDIR=`pwd`

export FDSUG=$SVNROOT/FDS/Manuals/FDS_User_Guide
export SMVUG=$SVNROOT/SMV/Manuals/SMV_User_Guide
export SMVVG=$SVNROOT/SMV/Manuals/SMV_Verification_Guide
SUMMARY=$SVNROOT/SMV/Manuals/SMV_Summary

is_file_installed $SMV
is_file_installed $SMOKEZIP
is_file_installed $SMOKEDIFF
is_file_installed $BACKGROUND
is_file_installed $WIND2FDS

cd $SMVUG/SCRIPT_FIGURES
rm -f *.png
rm -f *.help

rm -f smokeview.version
rm -f smokediff.version
rm -f smokezip.version
rm -f background.version
rm -f wind2fds.version

rm -f smokeview.help
rm -f smokediff.help
rm -f smokezip.help
rm -f background.help
rm -f wind2fds.help

rm -f $SUMMARY/images/*.png
source ~/.bashrc_fds

$SMV -help > smokeview.help
$SMOKEZIP -help > smokezip.help
$SMOKEDIFF -help > smokediff.help
$BACKGROUND -help > background.help
$WIND2FDS -help > wind2fds.help

$SMV -version > smokeview.version
$SMOKEZIP -v > smokezip.version
$SMOKEDIFF -v > smokediff.version
$BACKGROUND -version > background.version
$WIND2FDS  > wind2fds.version

cd $SMVVG/SCRIPT_FIGURES
rm -f *.version
rm -f *.png
$SMV -version > smokeview.version

if [ "$RUN_SMV" == "1" ] ; then
  cd $SVNROOT/SMV/Verification/Visualization
  echo Converting particles to isosurfaces in case plumeiso
  $SMOKEZIP -r -part2iso plumeiso

  cd $SVNROOT/SMV/Verification/WUI
  echo Converting particles to isosurfaces in case pine_tree
  if  [ -e pine_tree.smv ]; then
    $SMOKEZIP -r -part2iso pine_tree
  fi

# precompute FED slices

  source $STARTX 2>/dev/null
  $QFDS -f -d Visualization plume5c
  $QFDS -f -d Visualization plume5cdelta
  $QFDS -f -d Visualization thouse5
  $QFDS -f -d Visualization thouse5delta
  source $STOPX 2>/dev/null

# difference plume5c and thouse5

  cd $SVNROOT/SMV/Verification/Visualization
  echo Differencing cases plume5c and plume5cdelta
  $SMOKEDIFF -w -r plume5c plume5cdelta
  echo Differencing cases thouse5 and thouse5delta
  $SMOKEDIFF -w -r thouse5 thouse5delta

  echo Generating images

  source $STARTX
  cd $SVNROOT/SMV/Verification
  scripts/SMV_Cases.sh
  cd $SVNROOT/SMV/Verification
  scripts/SMV_DIFF_Cases.sh
  cd $CURDIDR
  source $STOPX

fi

# generate geometry images

if [ "$RUN_WUI" == "1" ] ; then
  source $STARTX
  cd $SVNROOT/SMV/Verification
  scripts/WUI_Cases.sh
  source $STOPX
fi
if [ "$RUN_GEOM" == "1" ] ; then
  source $STARTX
  cd $SVNROOT/SMV/Verification
  scripts/GEOM_Cases.sh
  source $STOPX
fi

# copy generated images to web summary directory

cp $SMVVG/FIGURES/graysquares.png $SUMMARY/images/.
cp $FDSUG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
cp $SMVUG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
cp $SMVVG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
