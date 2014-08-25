#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR PRISME_LK_1_Lower.fds
$QFDS  $QUEUE -d $INDIR PRISME_LK_1_Upper.fds
$QFDS  $QUEUE -d $INDIR PRISME_LK_2_Lower.fds
$QFDS  $QUEUE -d $INDIR PRISME_LK_2_Upper.fds
$QFDS  $QUEUE -d $INDIR PRISME_LK_3.fds      
$QFDS  $QUEUE -d $INDIR PRISME_LK_4.fds      

echo FDS cases submitted
