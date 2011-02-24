#!/bin/bash -f

export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/INTEL_LINUX_64/smokeview_linux_64
export SMOKEZIP=$SVNROOT/Utilities/smokezip/INTEL_LINUX_64/smokezip_linux_64
export SMOKEDIFF=$SVNROOT/Utilities/smokediff/INTEL_LINUX_64/smokediff_linux_64
export BACKGROUND=$SVNROOT/Utilities/background/INTEL_LINUX_32/background
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export BASEDIR=`pwd`

export SMVUG=$SVNROOT/Manuals/SMV_User_Guide
export SMVVG=$SVNROOT/Manuals/SMV_Verification_Guide


cd $SMVUG
rm SCRIPT_FIGURES/*.png
rm SCRIPT_FIGURES/*.help
rm SCRIPT_FIGURES/*.version

$SMV -help > SCRIPT_FIGURES\smokeview.help
$SMV -version > SCRIPT_FIGURES\smokeview.version
$SMOKEZIP -help > SCRIPT_FIGURES\smokezip.help
$SMOKEDIFF -help > SCRIPT_FIGURES\smokediff.help
$SMOKEDIFF -v > SCRIPT_FIGURES\smokediff.version
$BACKGROUND -help > SCRIPT_FIGURES\background.help
$BACKGROUND -version > SCRIPT_FIGURES\background.version

cd $SMVVG
rm SCRIPT_FIGURES\*.version
rm SCRIPT_FIGURES\*.png
$SMV -version > SCRIPT_FIGURES\smokeview.version

cd $SVNROOT/Verification
./SMV_Cases.sh
