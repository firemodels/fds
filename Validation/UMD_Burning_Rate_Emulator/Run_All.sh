#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_1_2D.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_2_2D.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_3_2D.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_4_2D.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_5_2D.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE_Test_6_2D.fds

echo FDS cases submitted
