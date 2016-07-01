#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_3.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_5.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_6.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_7.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_8.fds 
$QFDS $DEBUG $QUEUE -d $INDIR NRL_HAI_9.fds 

echo FDS cases submitted
