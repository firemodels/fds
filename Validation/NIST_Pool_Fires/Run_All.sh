#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 12 -n 8 $QUEUE -d $INDIR NIST_Acetone_Prescribed_2cm.fds
$QFDS $DEBUG -p 12 -n 8 $QUEUE -d $INDIR NIST_Acetone_Prescribed_1cm.fds
$QFDS $DEBUG -p 24 -n 8 $QUEUE -d $INDIR NIST_Acetone_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 12 -n 8 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 12 -n 8 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 24 -n 8 $QUEUE -d $INDIR NIST_Ethanol_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 16 -n 8 $QUEUE -d $INDIR NIST_Methanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 16 -n 8 $QUEUE -d $INDIR NIST_Methanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 32 -n 8 $QUEUE -d $INDIR NIST_Methanol_Prescribed_0p5cm.fds

$QFDS $DEBUG -p 96 -n 8 $QUEUE -d $INDIR NIST_Methanol_1m_pan_1cm_grid.fds
$QFDS $DEBUG -p 12 -n 6 $QUEUE -d $INDIR NIST_Methanol_1m_pan_2cm_grid.fds
$QFDS $DEBUG -p  2 -n 2 $QUEUE -d $INDIR NIST_Methanol_1m_pan_4cm_grid.fds
$QFDS $DEBUG -p 96 -n 4 $QUEUE -d $INDIR NIST_Methanol_1m_pan_1cm_grid_predicted.fds
$QFDS $DEBUG -p 12 -n 4 $QUEUE -d $INDIR NIST_Methanol_1m_pan_2cm_grid_predicted.fds
$QFDS $DEBUG -p  2 -n 2 $QUEUE -d $INDIR NIST_Methanol_1m_pan_4cm_grid_predicted.fds

echo FDS cases submitted
