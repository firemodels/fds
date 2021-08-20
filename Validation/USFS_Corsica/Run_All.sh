#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_1_0.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_2_0.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_3_0.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_1_20.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_2_20.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Test_3_20.fds

echo FDS cases submitted
