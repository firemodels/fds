#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 8 $DEBUG $QUEUE -d $INDIR tree_2_m_14_pc.fds
$QFDS -p 8 $DEBUG $QUEUE -d $INDIR tree_2_m_49_pc.fds
$QFDS -p 8 $DEBUG $QUEUE -d $INDIR tree_5_m_26_pc.fds

echo FDS cases submitted
