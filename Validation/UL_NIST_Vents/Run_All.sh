#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR UL_NIST_Vents_Test_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR UL_NIST_Vents_Test_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR UL_NIST_Vents_Test_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR UL_NIST_Vents_Test_4.fds

echo FDS cases submitted
