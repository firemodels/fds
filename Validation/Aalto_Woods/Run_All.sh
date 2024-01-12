#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR pine_flaming_25.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_flaming_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_flaming_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_N2_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_smolder_20.fds
$QFDS $DEBUG $QUEUE -d $INDIR pine_smolder_30.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_flaming_25.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_flaming_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_flaming_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_N2_35.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_N2_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_smolder_25.fds
$QFDS $DEBUG $QUEUE -d $INDIR spruce_smolder_35.fds
$QFDS -p 4 $DEBUG $QUEUE -d $INDIR Roomcorner_M12.fds
$QFDS -p 4 $DEBUG $QUEUE -d $INDIR Roomcorner_modified.fds

echo FDS cases submitted
