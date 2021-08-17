#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_post_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_post_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_post_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_post_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_pre_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_pre_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_pre_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Heil_pre_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_post_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_post_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_post_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_post_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_pre_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_pre_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_pre_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR Pike1_pre_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_post_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_post_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_post_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_post_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_pre_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_pre_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_pre_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR RedF_pre_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_post_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_post_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_post_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_post_9ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_pre_13ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_pre_2ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_pre_4ms.fds
$QFDS -p 72  $DEBUG $QUEUE -d $INDIR UNC_pre_9ms.fds

echo FDS cases submitted
