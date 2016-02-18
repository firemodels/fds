#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 8 $QUEUE -d $INDIR Smyth_Slot_Burner.fds
#$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Smyth_Slot_Burner_WD.fds
#$QFDS $DEBUG -p 32 $QUEUE -d $INDIR Smyth_Slot_Burner_Andersen.fds

echo FDS cases submitted
