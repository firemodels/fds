#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_C064_LS.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_F19_LS.fds

$QFDS -p 36 -n 6 $DEBUG $QUEUE -d $INDIR Case_C064.fds
$QFDS -p 36 -n 6 $DEBUG $QUEUE -d $INDIR Case_F19.fds

$QFDS -p 36 -n 6 $DEBUG $QUEUE -d $INDIR Case_C064_BFM.fds
$QFDS -p 36 -n 6 $DEBUG $QUEUE -d $INDIR Case_F19_BFM.fds

echo FDS cases submitted
