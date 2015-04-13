#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 12 $QUEUE -d $INDIR UMD_Gas_Burner_coarse.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR UMD_Gas_Burner_fine.fds

echo FDS cases submitted
