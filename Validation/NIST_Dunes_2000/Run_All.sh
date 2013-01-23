#!/bin/bash -f

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR NIST_Dunes_2000_SDC02
$RUNFDS $INDIR NIST_Dunes_2000_SDC05
$RUNFDS $INDIR NIST_Dunes_2000_SDC07
$RUNFDS $INDIR NIST_Dunes_2000_SDC10
$RUNFDS $INDIR NIST_Dunes_2000_SDC33
$RUNFDS $INDIR NIST_Dunes_2000_SDC35
$RUNFDS $INDIR NIST_Dunes_2000_SDC39
