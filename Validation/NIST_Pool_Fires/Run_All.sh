#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Acetone_Prescribed_2cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Acetone_Prescribed_1cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Acetone_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Methanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Methanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 12 $QUEUE -d $INDIR NIST_Methanol_Prescribed_0p5cm.fds

echo FDS cases submitted
