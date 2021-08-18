#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 65 $QUEUE -d $INDIR Askervein_TU03A_16m.fds
#$QFDS $DEBUG -p 145 $QUEUE -d $INDIR Askervein_TU03A_8m.fds
#$QFDS $DEBUG -p 145 $QUEUE -d $INDIR Askervein_TU03A_4m.fds

echo FDS cases submitted
