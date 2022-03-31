#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_3_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_4_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_8_4cm.fds
$QFDS $DEBUG $QUEUE -p 20 -d $INDIR Test_9_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_12_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_13_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_15_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_16_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_17_4cm.fds
$QFDS $DEBUG $QUEUE -p 40 -d $INDIR Test_19_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_22_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_23_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_24_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_26_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_27_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_29_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_30_4cm.fds
$QFDS $DEBUG $QUEUE -p 80 -d $INDIR Test_31_4cm.fds

