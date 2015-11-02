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
echo "-S host - make pictures on host"
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
SSH=

while getopts 'dghis:S:t' OPTION
do
case $OPTION  in
  d)
   DEBUG=_db
   ;;
  g)
   RUN_SMV=0
   RUN_GEOM=1
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
  S)
  SSH="ssh $OPTARG "
  ;;
esac
done
shift $(($OPTIND-1))


VERSION=$PLATFORM$TEST$SIZE$DEBUG
VERSION2=$PLATFORM$SIZE
CURDIR=`pwd`
cd ../..
export SVNROOT=`pwd`
cd $CURDIR/..

if [ "$use_installed" == "1" ] ; then
  export SMV=smokeview
  export SMOKEZIP=smokediff
  export SMOKEDIFF=smokediff
  export WIND2FDS=wind2fds
  export BACKGROUND=background
else
  export SMV=$SVNROOT/SMV/Build/intel_$VERSION2/smokeview_$VERSION
  export SMOKEZIP=$SVNROOT/Utilities/smokezip/intel_$VERSION2/smokezip_$VERSION2
  export SMOKEDIFF=$SVNROOT/Utilities/smokediff/intel_$VERSION2/smokediff_$VERSION2
  export WIND2FDS=$SVNROOT/Utilities/wind2fds/intel_$VERSION2/wind2fds_$VERSION2
  export BACKGROUND=$SVNROOT/Utilities/background/intel_$PLATFORM\_32/background
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
export RUNTFDS="$RUNSMV -t"
export RUNCFAST=$RUNSMV
export BASEDIR=`pwd`

export FDSUG=$SVNROOT/Manuals/FDS_User_Guide
export SMVUG=$SVNROOT/Manuals/SMV_User_Guide
export SMVVG=$SVNROOT/Manuals/SMV_Verification_Guide
SUMMARY=$SVNROOT/Manuals/SMV_Summary

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
  cd $SVNROOT/Verification/Visualization
  echo Converting particles to isosurfaces in case plumeiso
  $SMOKEZIP -r -part2iso plumeiso

  cd $SVNROOT/Verification/WUI
  echo Converting particles to isosurfaces in case pine_tree
  if  [ -e pine_tree.smv ]; then
    $SMOKEZIP -r -part2iso pine_tree
  fi

# precompute FED slices

  if [ "$SSH" == "" ]; then
  source $STARTX 2>/dev/null
  $QFDS -f -d Visualization plume5c
  $QFDS -f -d Visualization plume5cdelta
  $QFDS -f -d Visualization thouse5
  $QFDS -f -d Visualization thouse5delta
  source $STOPX 2>/dev/null
  else
  $SSH \( cd $SVNROOT/Verification \; source $STARTX 2>/dev/null \; \
  $QFDS -f -d Visualization plume5c \; \
  $QFDS -f -d Visualization plume5cdelta \; \ 
  $QFDS -f -d Visualization thouse5 \; \
  $QFDS -f -d Visualization thouse5delta \; \
  source $STOPX 2>/dev/null \)
  fi

# difference plume5c and thouse5

  cd $SVNROOT/Verification/Visualization
  echo Differencing cases plume5c and plume5cdelta
  $SMOKEDIFF -w -r plume5c plume5cdelta
  echo Differencing cases thouse5 and thouse5delta
  $SMOKEDIFF -w -r thouse5 thouse5delta

  echo Generating images

# copy wui error image in case wfds does not exist

  FROMDIR=$SVNROOT/Manuals/SMV_Verification_Guide/FIGURES
  TODIR=$SVNROOT/Manuals/SMV_Verification_Guide/SCRIPT_FIGURES
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_part_000.png
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_part_010.png
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_part_020.png
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_partiso_000.png
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_partiso_010.png
  cp $FROMDIR/wfds_error.png $TODIR/tree_one_partiso_020.png
 
  if [ "$SSH" == "" ]; then
  source $STARTX
  cd $SVNROOT/Verification
  scripts/SMV_Cases.sh
  cd $SVNROOT/Verification
  scripts/SMV_DIFF_Cases.sh
  cd $CURDIDR
  source $STOPX
  else
  $SSH \( source $STARTX \; \
  cd $SVNROOT/Verification \; \
  scripts/SMV_Cases.sh \; \
  cd $SVNROOT/Verification \; \
  scripts/SMV_DIFF_Cases.sh \; \
  cd $CURDIDR \; \
  source $STOPX \)
  fi

# copy generated images to web summary directory

  cp $SMVVG/FIGURES/graysquares.png $SUMMARY/images/.
  cp $FDSUG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
  cp $SMVUG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
  cp $SMVVG/SCRIPT_FIGURES/*.png $SUMMARY/images/.
fi

# generate geometry images

if [ "$RUN_GEOM" == "1" ] ; then
  if [ "$SSH" == "" ]; then
  source $STARTX
  cd $SVNROOT/Verification
  scripts/SMV_geom_Cases.sh
  source $STOPX
  else
  $SSH \( source $STARTX \; \
  cd $SVNROOT/Verification \; \
  scripts/SMV_geom_Cases.sh \; \
  source $STOPX \)
  fi
fi
