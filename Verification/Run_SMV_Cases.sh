#!/bin/bash -f

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

export SVNROOT=~/FDS-SMV
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export RUNFDSFG=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
# uncomment following line to stop all cases
#export STOPFDS=1

scripts/SMV_Cases.csh

echo FDS cases submitted

