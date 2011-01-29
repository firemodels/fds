#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR NRL_HAI_1
$RUNFDS $INDIR NRL_HAI_2
$RUNFDS $INDIR NRL_HAI_3 
$RUNFDS $INDIR NRL_HAI_4
$RUNFDS $INDIR NRL_HAI_5 
$RUNFDS $INDIR NRL_HAI_6 
$RUNFDS $INDIR NRL_HAI_7 
$RUNFDS $INDIR NRL_HAI_8 
$RUNFDS $INDIR NRL_HAI_9 

echo FDS cases submitted
