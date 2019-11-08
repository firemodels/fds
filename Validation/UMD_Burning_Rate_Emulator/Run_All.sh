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

$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_1_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_2_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_3_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_4_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_5_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_6_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_7_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_8_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_9_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_10_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_11_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_12_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_13_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_14_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_15_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_16_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_17_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_18_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_19_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_20_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_21_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_22_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_23_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_24_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_25_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_26_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_27_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_28_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_29_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_30_1mm.fds
$QFDS $DEBUG $QUEUE -d $INDIR BRE50_Test_31_1mm.fds

echo FDS cases submitted
