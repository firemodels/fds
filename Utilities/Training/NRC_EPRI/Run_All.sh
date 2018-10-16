#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Annulus_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Main_Control_Room_Purge_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Main_Control_Room_No_Purge_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR Pump_Room_v6.fds

$QFDS $DEBUG $QUEUE -d $INDIR -p 7 Cable_Spreading_Room_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 4 Corridor_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 5 Switchgear_Room_Cabinet_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 4 Switchgear_Room_MCC_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 2 Turbine_Building_Location_1_v6.fds
$QFDS $DEBUG $QUEUE -d $INDIR -p 2 Turbine_Building_Location_2_v6.fds

