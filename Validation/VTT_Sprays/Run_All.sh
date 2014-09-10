#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR LN02_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR LN02_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR LN02_4.fds

echo FDS cases submitted
