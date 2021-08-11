#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 160 $QUEUE -d $INDIR Montoir01.fds
$QFDS $DEBUG -p 160 $QUEUE -d $INDIR Montoir02.fds
$QFDS $DEBUG -p 160 $QUEUE -d $INDIR Montoir03.fds

echo FDS cases submitted
