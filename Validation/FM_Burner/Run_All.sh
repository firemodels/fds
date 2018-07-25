#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 125 -n 12 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_CH4.fds
$QFDS -p 125 -n 12 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C2H4.fds
$QFDS -p 125 -n 12 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C3H6.fds
$QFDS -p 125 -n 12 $DEBUG $QUEUE -d $INDIR FM_15cm_Burner_C3H8.fds

echo FDS cases submitted
