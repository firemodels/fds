#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Bakken_2p1.fds
$QFDS $DEBUG $QUEUE -p 9  -d $INDIR Tight_1.fds

echo FDS cases submitted
