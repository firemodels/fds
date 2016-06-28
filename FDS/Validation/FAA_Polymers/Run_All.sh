#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_HDPE.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_HIPS.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_35_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_50_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBT_70_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_35_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_50_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PBTGF_70_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PMMA.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PC_3_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PC_6_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PC_6_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PC_6_92.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PC_9_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PVC_3_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PVC_6_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PVC_6_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PVC_6_92.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PVC_9_75.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_90.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_50_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_70_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PEEK_90_solid_only.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PP.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PA66.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_POM.fds
$QFDS $DEBUG $QUEUE -d $INDIR FAA_Polymers_PET.fds
$QFDS $DEBUG $QUEUE -d $INDIR Cardboard_DW_20.fds
$QFDS $DEBUG $QUEUE -d $INDIR Cardboard_DW_40.fds
$QFDS $DEBUG $QUEUE -d $INDIR Cardboard_DW_60.fds

echo FDS cases submitted
