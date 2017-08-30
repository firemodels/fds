#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh


$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_1.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_2.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_3.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_4.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_5.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_6.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_7.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_8.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_9.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_10.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_11.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_12.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_13.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_14.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_15.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_16.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_17.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_18.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_19.fds
$QFDS $DEBUG $QUEUE -p 3 -n 3 -d $INDIR NFSC2_Room_20.fds

