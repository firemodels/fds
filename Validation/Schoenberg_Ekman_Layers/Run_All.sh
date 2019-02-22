#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR ekman_neutral.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR ekman_stable.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR ekman_unstable.fds

echo FDS cases submitted
