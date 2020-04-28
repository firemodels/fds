#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 96 -n 8 $QUEUE -d $INDIR PBT_60_kW.fds
$QFDS $DEBUG -p 96 -n 8 $QUEUE -d $INDIR PMMA_60_kW.fds
$QFDS $DEBUG -p 96 -n 8 $QUEUE -d $INDIR PVC_60_kW.fds

echo FDS cases submitted
