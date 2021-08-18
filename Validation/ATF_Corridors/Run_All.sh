#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_050_kW.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_100_kW.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_240_kW.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_250_kW.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_500_kW.fds
$QFDS $DEBUG $QUEUE -p 3 -d $INDIR ATF_Corridors_Mix_kW.fds
