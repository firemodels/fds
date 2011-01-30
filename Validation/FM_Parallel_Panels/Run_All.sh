#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR FM_Parallel_Panel_1 
$RUNFDS $INDIR FM_Parallel_Panel_2 
$RUNFDS $INDIR FM_Parallel_Panel_3 
$RUNFDS $INDIR FM_Parallel_Panel_4 
$RUNFDS $INDIR FM_Parallel_Panel_5 
$RUNFDS $INDIR FM_Parallel_Panel_6 

echo FDS cases submitted
