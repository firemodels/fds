#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR Hasemi_1D_9.fds
$QFDS -p 18 $DEBUG $QUEUE -d $INDIR Hasemi_1D_21.fds
$QFDS -p 54 $DEBUG $QUEUE -d $INDIR Hasemi_1D_33.fds

echo FDS cases submitted
