#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 12  $DEBUG $QUEUE -d $INDIR BST_FRS_1.fds
$QFDS -p 12  $DEBUG $QUEUE -d $INDIR BST_FRS_2.fds
$QFDS -p 12  $DEBUG $QUEUE -d $INDIR BST_FRS_3.fds
$QFDS -p 12  $DEBUG $QUEUE -d $INDIR BST_FRS_4.fds
$QFDS -p 12  $DEBUG $QUEUE -d $INDIR BST_FRS_5.fds
$QFDS -p 16  $DEBUG $QUEUE -d $INDIR BST_FRS_6.fds

echo FDS cases submitted
