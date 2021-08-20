#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 12 $DEBUG $QUEUE -d $INDIR CSTB_Tunnel_Test_2.fds
$QFDS -p 12 $DEBUG $QUEUE -d $INDIR CSTB_Tunnel_Test_27.fds

echo FDS cases submitted
