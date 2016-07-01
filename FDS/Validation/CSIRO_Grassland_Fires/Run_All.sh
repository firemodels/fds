#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 36 -n 3 $DEBUG $QUEUE -d $INDIR Case_C064.fds
$QFDS -p 36 -n 3 $DEBUG $QUEUE -d $INDIR Case_F19.fds

echo FDS cases submitted
