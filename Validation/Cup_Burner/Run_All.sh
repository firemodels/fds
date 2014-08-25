#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  -p 4 $QUEUE -d $INDIR Cup_C7H16_Ar.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_C7H16_CO2.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_C7H16_He.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_C7H16_N2.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_CH4_Ar.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_CH4_CO2.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_CH4_He.fds
$QFDS  -p 4 $QUEUE -d $INDIR Cup_CH4_N2.fds

echo FDS cases submitted
