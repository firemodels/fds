#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_21.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_22.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_23.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_24.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_42.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_46.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_55.fds
$QFDS $DEBUG -p 8 -r $QUEUE -d $INDIR Prairie_Grass_Case_57.fds 

echo FDS cases submitted
