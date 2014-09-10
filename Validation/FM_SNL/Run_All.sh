#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_01.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_02.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_03.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_04.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_05.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_06.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_07.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_08.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_09.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_10.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_11.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_12.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_13.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_14.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_15.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR FM_SNL_16.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR FM_SNL_17.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_21.fds
$QFDS $DEBUG -p 5 $QUEUE -d $INDIR FM_SNL_22.fds
