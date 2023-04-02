#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R3.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R4.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Gasification_Test_R5.fds

echo FDS cases submitted
