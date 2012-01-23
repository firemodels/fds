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
$RUNFDS $INDIR FAA_Polymers_PBT
$RUNFDS $INDIR FAA_Polymers_PBT_35
$RUNFDS $INDIR FAA_Polymers_PBT_50
$RUNFDS $INDIR FAA_Polymers_PBT_70
$RUNFDS $INDIR FAA_Polymers_PBT_35_solid_only
$RUNFDS $INDIR FAA_Polymers_PBT_50_solid_only
$RUNFDS $INDIR FAA_Polymers_PBT_70_solid_only
$RUNFDS $INDIR FAA_Polymers_PBTGF
$RUNFDS $INDIR FAA_Polymers_PBTGF_35
$RUNFDS $INDIR FAA_Polymers_PBTGF_50
$RUNFDS $INDIR FAA_Polymers_PBTGF_70
$RUNFDS $INDIR FAA_Polymers_PBTGF_35_solid_only
$RUNFDS $INDIR FAA_Polymers_PBTGF_50_solid_only
$RUNFDS $INDIR FAA_Polymers_PBTGF_70_solid_only
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
$RUNFDS $INDIR FAA_Polymers_PEEK_50
$RUNFDS $INDIR FAA_Polymers_PEEK_70
$RUNFDS $INDIR FAA_Polymers_PEEK_90
$RUNFDS $INDIR FAA_Polymers_PEEK_50_solid_only
$RUNFDS $INDIR FAA_Polymers_PEEK_70_solid_only
$RUNFDS $INDIR FAA_Polymers_PEEK_90_solid_only

echo FDS cases submitted
