#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 24 -d $INDIR BK-R.fds
$QFDS $DEBUG $QUEUE -p 24 -d $INDIR CLC-I-R.fds
$QFDS $DEBUG $QUEUE -p 24 -d $INDIR CLC-II-R.fds
$QFDS $DEBUG $QUEUE -p 24 -d $INDIR CMP-R.fds
