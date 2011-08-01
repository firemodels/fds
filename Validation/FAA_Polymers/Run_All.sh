#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR FAA_Polymers_HDPE
$RUNFDS $INDIR FAA_Polymers_HIPS
$RUNFDS $INDIR FAA_Polymers_PMMA
$RUNFDS $INDIR FAA_Polymers_PC_3_75
$RUNFDS $INDIR FAA_Polymers_PC_6_50
$RUNFDS $INDIR FAA_Polymers_PC_6_75
$RUNFDS $INDIR FAA_Polymers_PC_6_92
$RUNFDS $INDIR FAA_Polymers_PC_9_75
$RUNFDS $INDIR FAA_Polymers_PVC_3_75
$RUNFDS $INDIR FAA_Polymers_PVC_6_50
$RUNFDS $INDIR FAA_Polymers_PVC_6_75
$RUNFDS $INDIR FAA_Polymers_PVC_6_92
$RUNFDS $INDIR FAA_Polymers_PVC_9_75

echo FDS cases submitted
