#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_01.fds
$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_05.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_07.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_19.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_21.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Hamins_CH4_23.fds 

echo FDS cases submitted
