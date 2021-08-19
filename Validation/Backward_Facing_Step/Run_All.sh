#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG       $QUEUE -d $INDIR backward_facing_step_5.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR backward_facing_step_10.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR backward_facing_step_20.fds

echo FDS cases submitted
