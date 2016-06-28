#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_2.fds 
$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_3.fds 
$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_4.fds 
$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_5.fds 
$QFDS $DEBUG $QUEUE -d $INDIR FM_Parallel_Panel_6.fds 

echo FDS cases submitted
