#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Bryant_034_kW
$RUNFDS $INDIR Bryant_065_kW
$RUNFDS $INDIR Bryant_096_kW
$RUNFDS $INDIR Bryant_128_kW
$RUNFDS $INDIR Bryant_160_kW
$RUNFDS $INDIR Bryant_320_kW
$RUNFDS $INDIR Bryant_511_kW
 
echo FDS cases submitted
