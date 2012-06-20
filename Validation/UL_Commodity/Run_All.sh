#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Group_A_2x2x2
$RUNFDS $INDIR Group_A_2x2x3 
$RUNFDS $INDIR Group_A_2x2x4
$RUNFDS $INDIR Group_A_FM_RDD_p21.fds
$RUNFDS $INDIR Group_A_FM_RDD_p31.fds
$RUNFDS $INDIR Group_A_FM_RDD_p39.fds

echo FDS cases submitted
