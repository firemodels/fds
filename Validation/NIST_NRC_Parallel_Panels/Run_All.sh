#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 96  $QUEUE -d $INDIR PBT_60_kW.fds
$QFDS $DEBUG -p 12  $QUEUE -d $INDIR PMMA_60_kW_4_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PMMA_60_kW_2_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR PMMA_60_kW_1_cm.fds
$QFDS $DEBUG -p 96  $QUEUE -d $INDIR PVC_60_kW.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Marinite_60_kW_2_cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Marinite_60_kW_1_cm.fds
$QFDS $DEBUG -p 288 $QUEUE -d $INDIR Marinite_60_kW_5_mm.fds

echo FDS cases submitted
