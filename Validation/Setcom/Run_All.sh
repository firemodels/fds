#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 3 $QUEUE -d $INDIR setcom_0p8.fds
$QFDS $DEBUG -p 3 $QUEUE -d $INDIR setcom_1p8.fds
$QFDS $DEBUG -p 3 $QUEUE -d $INDIR setcom_3p7.fds
$QFDS $DEBUG -p 3 $QUEUE -d $INDIR setcom_4p2.fds
$QFDS $DEBUG -p 3 $QUEUE -d $INDIR setcom_5p2.fds

echo FDS cases submitted
