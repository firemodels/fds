#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 5 -n 5 $QUEUE -d $INDIR  FAA_B707_Test_1.fds
$QFDS $DEBUG -p 5 -n 5 $QUEUE -d $INDIR  FAA_B707_Test_2.fds
$QFDS $DEBUG -p 5 -n 5 $QUEUE -d $INDIR  FAA_B707_Test_3.fds

$QFDS $DEBUG -p  12 -n 6 $QUEUE -d $INDIR  FAA_B747_geom_6cm.fds
$QFDS $DEBUG -p 144 -n 6 $QUEUE -d $INDIR  FAA_B747_geom_4cm.fds
#$QFDS $DEBUG -p 144 -n 6 $QUEUE -d $INDIR  FAA_B747_geom_2cm.fds
