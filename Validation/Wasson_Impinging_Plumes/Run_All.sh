#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 75  $QUEUE -d $INDIR wasson_test1_12mm.fds
$QFDS $DEBUG -p 75  $QUEUE -d $INDIR wasson_test2_12mm.fds
$QFDS $DEBUG -p 50  $QUEUE -d $INDIR wasson_test3_12mm.fds
$QFDS $DEBUG -p 50  $QUEUE -d $INDIR wasson_test3_25mm.fds
$QFDS $DEBUG -p 9   $QUEUE -d $INDIR wasson_test3_50mm.fds
$QFDS $DEBUG -p 100 $QUEUE -d $INDIR wasson_test4_12mm.fds
$QFDS $DEBUG -p 75  $QUEUE -d $INDIR wasson_test5_12mm.fds
$QFDS $DEBUG -p 75  $QUEUE -d $INDIR wasson_test6_12mm.fds

echo FDS cases submitted
