#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_20kg_050cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_20kg_100cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_20kg_150cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_20kg_200cms.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_40kg_050cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_40kg_100cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_40kg_150cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_40kg_200cms.fds

$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_60kg_050cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_60kg_100cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_60kg_150cms.fds
$QFDS $DEBUG -p 4 $QUEUE -d $INDIR PineNeedles_60kg_200cms.fds

echo FDS cases submitted
