#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 12  $QUEUE -d $INDIR PBT_4_cm.fds
$QFDS $DEBUG -p 72  $QUEUE -d $INDIR PBT_2_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PBT_1_cm.fds
$QFDS $DEBUG -p 12  $QUEUE -d $INDIR PMMA_60_kW_4_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PMMA_60_kW_2_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PMMA_60_kW_1_cm.fds
$QFDS $DEBUG -p 12  $QUEUE -d $INDIR PVC_4_cm.fds
$QFDS $DEBUG -p 72  $QUEUE -d $INDIR PVC_2_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PVC_1_cm.fds
$QFDS $DEBUG -p 8   $QUEUE -d $INDIR Marinite_60_kW_4_cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Marinite_60_kW_2_cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Marinite_60_kW_1_cm.fds
$QFDS $DEBUG -p 12  $QUEUE -d $INDIR WRC_4_cm.fds
$QFDS $DEBUG -p 72  $QUEUE -d $INDIR WRC_2_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR WRC_1_cm.fds

echo FDS cases submitted
