#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Hamins_CH4_01 
$RUNFDS $INDIR Hamins_CH4_05 
$RUNFDS $INDIR Hamins_CH4_07 
$RUNFDS $INDIR Hamins_CH4_19 
$RUNFDS $INDIR Hamins_CH4_21 
$RUNFDS $INDIR Hamins_CH4_23 

echo FDS cases submitted
