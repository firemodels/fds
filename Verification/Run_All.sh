#!/bin/bash -f
export SVNROOT=~/FDS-SMV
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds.sh
export SMOKEZIP=$SVNROOT/Utilities/smokezip/INTEL_LINUX_32/smokezip_linux_32
export BASEDIR=`pwd`

./FDS_Cases.csh

echo FDS cases submitted

cd $BASEDIR/scripts

./SMV_Cases.csh

echo SMV cases submitted

cd $BASEDIR/scripts

./run_wui_tree_test.csh

echo WUI case submitted
