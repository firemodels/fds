#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Waterloo_Methanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Waterloo_Methanol_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Waterloo_Methanol_Prescribed_0p25cm.fds
#$QFDS $DEBUG -p 128 -n 8 $QUEUE -d $INDIR Waterloo_Methanol_Prescribed_0p125cm.fds

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR Waterloo_Methanol_Predicted_1cm.fds
$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Waterloo_Methanol_Predicted_0p5cm.fds
$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Waterloo_Methanol_Predicted_0p25cm.fds
#$QFDS $DEBUG -p 128 -n 8 $QUEUE -d $INDIR Waterloo_Methanol_Predicted_0p125cm.fds

echo FDS cases submitted
