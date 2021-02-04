#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 30 $DEBUG $QUEUE -d $INDIR Wu_Bakar_Tunnel_A.fds
$QFDS -p 30 $DEBUG $QUEUE -d $INDIR Wu_Bakar_Tunnel_B.fds
$QFDS -p 30 $DEBUG $QUEUE -d $INDIR Wu_Bakar_Tunnel_C.fds
$QFDS -p 30 $DEBUG $QUEUE -d $INDIR Wu_Bakar_Tunnel_D.fds
$QFDS -p 30 $DEBUG $QUEUE -d $INDIR Wu_Bakar_Tunnel_E.fds

echo FDS cases submitted
