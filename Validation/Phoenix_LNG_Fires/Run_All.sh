#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 34 $QUEUE -d $INDIR Phoenix01.fds
$QFDS $DEBUG -p 18 $QUEUE -d $INDIR Phoenix02.fds

echo FDS cases submitted
