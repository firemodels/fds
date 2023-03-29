#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 36  $QUEUE -d $INDIR UMD_SBI_1_cm.fds
$QFDS $DEBUG -p 128 $QUEUE -d $INDIR UMD_SBI_5_mm.fds

echo FDS cases submitted
