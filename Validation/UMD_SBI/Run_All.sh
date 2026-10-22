#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 8  $QUEUE -d $INDIR UMD_SBI_4_cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR UMD_SBI_2_cm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR UMD_SBI_1_cm.fds

echo FDS cases submitted
