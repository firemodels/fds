#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR ABS_30.fds
$QFDS $DEBUG $QUEUE -d $INDIR ABS_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR ABS_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR HIPS_30.fds
$QFDS $DEBUG $QUEUE -d $INDIR HIPS_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR HIPS_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR Kydex_30.fds
$QFDS $DEBUG $QUEUE -d $INDIR Kydex_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR Kydex_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR PEI_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR PEI_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR PEI_90.fds
$QFDS $DEBUG $QUEUE -d $INDIR PET_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR PET_70.fds
$QFDS $DEBUG $QUEUE -d $INDIR PMMA_20.fds
$QFDS $DEBUG $QUEUE -d $INDIR PMMA_40.fds
$QFDS $DEBUG $QUEUE -d $INDIR PMMA_60.fds
$QFDS $DEBUG $QUEUE -d $INDIR POM_30.fds
$QFDS $DEBUG $QUEUE -d $INDIR POM_50.fds
$QFDS $DEBUG $QUEUE -d $INDIR POM_70.fds

echo FDS cases submitted
