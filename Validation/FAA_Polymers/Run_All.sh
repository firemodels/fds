#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR FAA_Polymers_HDPE.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_HIPS.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_35.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_50.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_70.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_35_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_50_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBT_70_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_35.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_50.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_70.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_35_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_50_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PBTGF_70_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PMMA.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PC_3_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PC_6_50.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PC_6_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PC_6_92.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PC_9_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PVC_3_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PVC_6_50.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PVC_6_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PVC_6_92.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PVC_9_75.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_50.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_70.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_90.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_50_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_70_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PEEK_90_solid_only.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PP.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PA66.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_POM.fds
$QFDS -r $qq -d $INDIR FAA_Polymers_PET.fds

echo FDS cases submitted
