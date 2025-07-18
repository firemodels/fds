#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 15 $QUEUE -d $INDIR FM_Datacenter_Veltest_Low.fds
$QFDS $DEBUG -p 15 $QUEUE -d $INDIR FM_Datacenter_Veltest_High.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_Low_C3H6_SF.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_Low_C3H6_HA.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_High_C3H6_SF.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_High_C3H6_HA.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_Low_Cable_SF.fds
$QFDS $DEBUG -p 17 $QUEUE -d $INDIR FM_Datacenter_High_Cable_SF.fds

echo FDS cases submitted
