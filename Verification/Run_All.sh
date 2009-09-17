#!/bin/bash -f
export SVNROOT=~/FDS-SMV
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds5_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`

./FDS_Cases.csh

echo FDS cases submitted
