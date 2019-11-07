#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 10 -n 10 $DEBUG $QUEUE -d $INDIR SP_Test_3.fds
$QFDS -p 10 -n 10 $DEBUG $QUEUE -d $INDIR SP_Test_6.fds
$QFDS -p 10 -n 10 $DEBUG $QUEUE -d $INDIR SP_Test_9.fds

echo FDS cases submitted
