#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR Ulster_SBI_30_kW.fds
$QFDS  $QUEUE -d $INDIR Ulster_SBI_45_kW.fds
$QFDS  $QUEUE -d $INDIR Ulster_SBI_60_kW.fds 

echo FDS cases submitted
