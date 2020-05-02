#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

cp -v ./FDS_Input_Files_GEOM/*.fds $INDIR

$QFDS $DEBUG -p  21 -n 12 $QUEUE -d $INDIR  FAA_B747_front11kw_6cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_front11kw_4cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_front11kw_2cm.fds
#$QFDS $DEBUG -p 288 -n 12 $QUEUE -d $INDIR  FAA_B747_front11kw_1p25cm.fds
$QFDS $DEBUG -p  21 -n 12 $QUEUE -d $INDIR  FAA_B747_front5p5kw_6cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_front5p5kw_4cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_front5p5kw_2cm.fds
$QFDS $DEBUG -p  21 -n 12 $QUEUE -d $INDIR  FAA_B747_rear11kw_6cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_rear11kw_4cm.fds
$QFDS $DEBUG -p 234 -n 12 $QUEUE -d $INDIR  FAA_B747_rear11kw_2cm.fds

