#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 1 -o 8 -d $INDIR hu_1p.fds
$QFDS $DEBUG $QUEUE -p 1 -o 8 -d $INDIR hu_2p.fds
$QFDS $DEBUG $QUEUE -p 1 -o 8 -d $INDIR hu_3p.fds
$QFDS $DEBUG $QUEUE -p 1 -o 8 -d $INDIR hu_4p.fds

echo FDS cases submitted
