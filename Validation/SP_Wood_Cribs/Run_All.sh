#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_1.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_3.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_4.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_5.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_6.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_7.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_8.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_9.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_10.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_11.fds
$QFDS -p 10  $DEBUG $QUEUE -d $INDIR SP_Test_12.fds

echo FDS cases submitted
