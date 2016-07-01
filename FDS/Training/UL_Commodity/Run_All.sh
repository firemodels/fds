#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Group_A_2x2x2.fds
$QFDS $DEBUG $QUEUE -d $INDIR Group_A_2x2x3.fds
$QFDS $DEBUG $QUEUE -d $INDIR Group_A_2x2x4.fds
$QFDS $DEBUG $QUEUE -d $INDIR Group_A_FM_RDD_p21.fds
$QFDS $DEBUG $QUEUE -d $INDIR Group_A_FM_RDD_p31.fds
$QFDS $DEBUG $QUEUE -d $INDIR Group_A_FM_RDD_p39.fds

echo FDS cases submitted
