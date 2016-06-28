#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR neutral_bl.fds 
$QFDS $DEBUG $QUEUE -d $INDIR conv_bl.fds 
$QFDS $DEBUG $QUEUE -d $INDIR stable_bl.fds 

echo FDS cases submitted
