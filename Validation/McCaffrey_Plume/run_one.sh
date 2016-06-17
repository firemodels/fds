#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_Forty_Fives.sh

$QFDS $DEBUG -p 27 $QUEUE -d $INDIR McCaffrey_57_kW_45.fds

echo FDS cases submitted
