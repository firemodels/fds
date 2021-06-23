#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF1_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF2_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF3_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF4_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF5_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF6_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF7_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF8_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF9_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF10_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF11_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF12_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF13_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF14_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF15_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF16_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF17_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF18_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF19_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF20_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF21_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF22_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF25_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF26_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF27_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF28_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF29_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF30_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF31_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF32_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF33_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF34_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF35_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF36_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF37_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF38_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF39_3cm.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR DSF40_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF41_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF42_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF43_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF44_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF45_3cm.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR DSF46_3cm.fds

echo FDS cases submitted
