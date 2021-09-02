#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 75 -d $INDIR Exp_2.fds
$QFDS $DEBUG $QUEUE -p 75 -d $INDIR Exp_3.fds
$QFDS $DEBUG $QUEUE -p 75 -d $INDIR Exp_4.fds

echo FDS cases submitted
