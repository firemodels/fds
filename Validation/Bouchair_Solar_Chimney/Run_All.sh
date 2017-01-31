#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_1p_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_1p_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_1p_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_1p_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p1_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p1_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p1_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p1_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p2_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p2_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p2_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p2_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p3_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p3_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p3_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p3_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p5_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p5_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p5_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p1_p5_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_1p_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_1p_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_1p_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_1p_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p1_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p1_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p1_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p1_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p2_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p2_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p2_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p2_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p3_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p3_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p3_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p3_60.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p5_30.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p5_40.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p5_50.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR SC_p4_p5_60.fds

