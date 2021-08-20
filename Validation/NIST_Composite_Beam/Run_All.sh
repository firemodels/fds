#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 32 -d $INDIR NIST_CB_Test_2.fds
$QFDS $DEBUG $QUEUE -p 32 -d $INDIR NIST_CB_Test_3.fds
$QFDS $DEBUG $QUEUE -p 32 -d $INDIR NIST_CB_Test_4.fds
$QFDS $DEBUG $QUEUE -p 32 -d $INDIR NIST_CB_Test_5.fds

echo FDS cases submitted
