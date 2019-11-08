#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_1_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_2_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_3_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_4_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_5_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_6_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_7_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_8_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_9_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_10_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_11_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_12_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_13_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_14_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_15_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_16_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_17_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE25_Test_18_1mm.fds

echo FDS cases submitted
