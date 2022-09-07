#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 48 $DEBUG $QUEUE -d $INDIR SSE_Verification_4MW.fds
$QFDS -p 48 $DEBUG $QUEUE -d $INDIR SSE_Verification_8MWb.fds

$QFDS -p 48 $DEBUG $QUEUE -d $INDIR SSE_Verification_4MW_10cm.fds
$QFDS -p 48 $DEBUG $QUEUE -d $INDIR SSE_Verification_8MWb_10cm.fds

echo FDS cases submitted
