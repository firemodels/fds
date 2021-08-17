#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_1.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_2.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_3.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_4.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_5.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_6.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_7.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_8.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_9.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_10.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_11.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_12.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_13.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_14.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_15.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_16.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_17.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_18.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_19.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_20.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_21.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_22.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_23.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_24.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_25.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_26.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_27.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_28.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_29.fds
$QFDS $DEBUG  $QUEUE -d $INDIR stone_wool_const_AEn_30.fds

echo FDS cases submitted
