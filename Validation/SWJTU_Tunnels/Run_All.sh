#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-1.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-2.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-3.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-4.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-5.fds
$QFDS $DEBUG $QUEUE -p 64 -d $INDIR Test_I-6.fds
