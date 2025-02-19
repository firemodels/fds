#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test1_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test1_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test2_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test2_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test3_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test3_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test4_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test4_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test5_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test5_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test6_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test6_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test7_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test7_jet_050mm.fds
$QFDS $DEBUG -p 14 $QUEUE -d $INDIR lattimer_test8_jet_025mm.fds
$QFDS $DEBUG -p 7  $QUEUE -d $INDIR lattimer_test8_jet_050mm.fds

echo FDS cases submitted


