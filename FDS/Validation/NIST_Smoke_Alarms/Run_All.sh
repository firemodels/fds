#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

# Run FDS cases
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC02.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC05.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC07.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC10.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC33.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC35.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC38.fds
$QFDS $DEBUG $QUEUE -d $INDIR NIST_Smoke_Alarms_SDC39.fds

echo FDS cases submitted
