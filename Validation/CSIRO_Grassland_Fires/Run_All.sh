#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_C064_LS.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_F19_LS.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_C064_LS_fine.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_F19_LS_fine.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_C064_LS_crude.fds
$QFDS -p 1 $DEBUG $QUEUE -d $INDIR Case_F19_LS_crude.fds

$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_F19.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064_fine.fds
$QFDS -p 144 $DEBUG $QUEUE -d $INDIR Case_F19_fine.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064_crude.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_F19_crude.fds

$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064_BFM.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_F19_BFM.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064_BFM_fine.fds
$QFDS -p 144 $DEBUG $QUEUE -d $INDIR Case_F19_BFM_fine.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_C064_BFM_crude.fds
$QFDS -p 36  $DEBUG $QUEUE -d $INDIR Case_F19_BFM_crude.fds

echo FDS cases submitted
