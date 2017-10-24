#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_1_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_2_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_3_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_4_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_5_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_6_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_7_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_8_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_9_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_10_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_11_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_12_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_13_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_14_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_15_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_16_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconv_17_8.fds

$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_1_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_2_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_3_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_4_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_5_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_6_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_7_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_8_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_9_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_10_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_11_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_12_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_13_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_14_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_15_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_16_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_17_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_18_8.fds
$QFDS $DEBUG -I $QUEUE -d $INDIR natconh_19_8.fds

$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_1_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_2_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_3_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_4_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_5_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_6_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_7_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_8_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_9_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_10_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_11_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_12_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_13_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_14_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_15_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_16_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconv_17_16.fds

$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_1_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_2_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_3_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_4_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_5_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_6_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_7_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_8_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_9_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_10_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_11_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_12_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_13_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_14_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_15_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_16_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_17_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_18_16.fds
$QFDS $DEBUG -I -p 16 -n 4 $QUEUE -d $INDIR natconh_19_16.fds

$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_1_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_2_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_3_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_4_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_5_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_6_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_7_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_8_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_9_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_10_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_11_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_12_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_13_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_14_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_15_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_16_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconv_17_32.fds

$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_1_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_2_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_3_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_4_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_5_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_6_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_7_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_8_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_9_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_10_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_11_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_12_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_13_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_14_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_15_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_16_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_17_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_18_32.fds
$QFDS $DEBUG -I -p 64 -n 4 $QUEUE -d $INDIR natconh_19_32.fds

echo FDS cases submitted
