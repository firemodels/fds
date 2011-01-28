#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR

$RUNFDS $INDIR McCaffrey_14_kW 
$RUNFDS $INDIR McCaffrey_22_kW 
$RUNFDS $INDIR McCaffrey_33_kW 
$RUNFDS $INDIR McCaffrey_45_kW 
$RUNFDS $INDIR McCaffrey_57_kW 

echo FDS cases submitted
