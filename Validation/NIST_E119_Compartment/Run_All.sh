#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 72 -d $INDIR E119_Compartment_Test_1.fds
$QFDS $DEBUG $QUEUE -p 72 -d $INDIR E119_Compartment_Test_2.fds
$QFDS $DEBUG $QUEUE -p 72 -d $INDIR E119_Compartment_Test_3.fds

echo FDS cases submitted
