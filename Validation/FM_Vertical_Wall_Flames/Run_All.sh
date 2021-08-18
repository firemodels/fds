#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 30 -d $INDIR propylene.fds
$QFDS $DEBUG $QUEUE -p 30 -d $INDIR ethane.fds
$QFDS $DEBUG $QUEUE -p 30 -d $INDIR ethylene.fds
$QFDS $DEBUG $QUEUE -p 30 -d $INDIR methane.fds

echo FDS cases submitted
