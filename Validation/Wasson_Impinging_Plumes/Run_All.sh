#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

#$QFDS $DEBUG -p 4  $QUEUE -d $INDIR wasson_50kW_1m_50mm.fds
#$QFDS $DEBUG -p 4  $QUEUE -d $INDIR wasson_50kW_1m_25mm.fds
$QFDS $DEBUG -p 32 $QUEUE -d $INDIR wasson_50kW_1m_12mm.fds

echo FDS cases submitted
