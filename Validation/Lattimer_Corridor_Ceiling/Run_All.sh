#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test1_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test1_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test1_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test2_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test2_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test2_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test3_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test3_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test3_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test4_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test4_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test4_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test5_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test5_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test5_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test6_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test6_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test6_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test7_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test7_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test7_050mm.fds
$QFDS $DEBUG -p 64 $QUEUE -d $INDIR lattimer_test8_012mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR lattimer_test8_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test8_050mm.fds

echo FDS cases submitted


