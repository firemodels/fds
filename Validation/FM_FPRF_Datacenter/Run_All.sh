#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Cable_Sub_H.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Cable_Sub_H_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Cable_Sub_L.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Cable_Sub_L_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_High.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_High_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Low.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Low_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Prop_Hot_H.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Prop_Hot_H_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Prop_Hot_L.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Prop_Hot_L_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Prop_Sub_H.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Prop_Sub_H_2.fds
$QFDS $DEBUG -r -p 2 $QUEUE -d $INDIR FM_Datacenter_Prop_Sub_L.fds
$QFDS $DEBUG -r -p 4 $QUEUE -d $INDIR FM_Datacenter_Prop_Sub_L_2.fds

echo FDS cases submitted
