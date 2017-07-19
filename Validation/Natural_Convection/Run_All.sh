#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR natconv_1_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_2_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_3_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_4_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_5_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_6_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_7_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_8_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_9_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_10_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_11_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_12_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconv_13_8.fds

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_1_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_2_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_3_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_4_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_5_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_6_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_7_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_8_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_9_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_10_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_11_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_12_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconv_13_16.fds

$QFDS $DEBUG $QUEUE -d $INDIR natconh_1_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_2_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_3_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_4_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_5_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_6_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_7_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_8_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_9_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_10_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_11_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_12_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_13_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_14_8.fds
$QFDS $DEBUG $QUEUE -d $INDIR natconh_15_8.fds

$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_1_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_2_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_3_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_4_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_5_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_6_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_7_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_8_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_9_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_10_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_11_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_12_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_13_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_14_16.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR natconh_15_16.fds

echo FDS cases submitted
