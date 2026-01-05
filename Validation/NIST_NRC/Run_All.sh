#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_01.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_02.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_03.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_04.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_05.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_07.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_08.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_09.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_10.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_13.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_14.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_15.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_16.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_17.fds
$QFDS $DEBUG -p 15  $QUEUE -d $INDIR NIST_NRC_18.fds

echo FDS cases submitted
