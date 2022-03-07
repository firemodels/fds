#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 18 $QUEUE -d $INDIR heated_channel_Pr_0p10_16.fds
$QFDS $DEBUG -p 18 $QUEUE -d $INDIR heated_channel_Pr_0p71_16.fds
$QFDS $DEBUG -p 18 $QUEUE -d $INDIR heated_channel_Pr_1p00_16.fds
$QFDS $DEBUG -p 18 $QUEUE -d $INDIR heated_channel_Pr_2p00_16.fds

echo FDS cases submitted
