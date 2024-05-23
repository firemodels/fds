#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR pine_10p5O2_25_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_10p5O2_25_3C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_10p5O2_40_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_10p5O2_40_3C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_21O2_25_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_21O2_25_3C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_21O2_40_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_21O2_40_3C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_25_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_25_3C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_40_1C.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_40_3C.fds

echo FDS cases submitted
