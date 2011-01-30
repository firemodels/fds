#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR WTC_01 
$RUNFDS $INDIR WTC_02 
$RUNFDS $INDIR WTC_03 
$RUNFDS $INDIR WTC_04 
$RUNFDS $INDIR WTC_05 
$RUNFDS $INDIR WTC_06 

echo FDS cases submitted
