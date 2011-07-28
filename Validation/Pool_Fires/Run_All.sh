#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR acetone_1_m
$RUNFDS $INDIR ethanol_1_m
$RUNFDS $INDIR methanol_1_m
$RUNFDS $INDIR butane_1_m
$RUNFDS $INDIR benzene_1_m
$RUNFDS $INDIR heptane_1_m


echo FDS cases submitted
