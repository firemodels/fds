#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR WTC_01.fds 
$QFDS $DEBUG $QUEUE -d $INDIR WTC_02.fds
$QFDS $DEBUG $QUEUE -d $INDIR WTC_03.fds 
$QFDS $DEBUG $QUEUE -d $INDIR WTC_04.fds 
$QFDS $DEBUG $QUEUE -d $INDIR WTC_05.fds 
$QFDS $DEBUG $QUEUE -d $INDIR WTC_06.fds 

echo FDS cases submitted
