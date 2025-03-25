#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 64  -d $INDIR crib_1x1x1_2cm.fds
$QFDS $DEBUG $QUEUE -p  8  -d $INDIR crib_1x1x1_4cm.fds
$QFDS $DEBUG $QUEUE -p  1  -d $INDIR crib_1x1x1_8cm.fds
$QFDS $DEBUG $QUEUE -p  4  -d $INDIR crib_2x1x1_8cm.fds
$QFDS $DEBUG $QUEUE -p  8  -d $INDIR crib_2x2x1_8cm.fds
$QFDS $DEBUG $QUEUE -p 12  -d $INDIR crib_2x2x2_8cm.fds

$QFDS $DEBUG $QUEUE -p 12  -d $INDIR box_1x1x1.fds
$QFDS $DEBUG $QUEUE -p 24  -d $INDIR box_2x1x1.fds
$QFDS $DEBUG $QUEUE -p 36  -d $INDIR box_2x2x1.fds
$QFDS $DEBUG $QUEUE -p 36  -d $INDIR box_2x2x2.fds

$QFDS $DEBUG $QUEUE -p 16  -d $INDIR pallet_1x1x2.fds
$QFDS $DEBUG $QUEUE -p 20  -d $INDIR pallet_1x1x4.fds
$QFDS $DEBUG $QUEUE -p 24  -d $INDIR pallet_1x1x8.fds

echo FDS cases submitted
