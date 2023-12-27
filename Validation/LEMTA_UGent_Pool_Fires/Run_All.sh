#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 16  $DEBUG $QUEUE -d $INDIR H5_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR H8_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR H10_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR H15_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR M5_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR M8_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR M10_1cm.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR M15_1cm.fds

echo FDS cases submitted
