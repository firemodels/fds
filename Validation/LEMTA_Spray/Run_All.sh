#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_3.fds
$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_4.fds
$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_5.fds
$QFDS $DEBUG $QUEUE -d $INDIR LEMTA_Spray_Cooling_SU42_Horizontal.fds

echo FDS cases submitted

