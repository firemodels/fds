#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Acetone_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Acetone_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Acetone_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Ethanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Ethanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Ethanol_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR NIST_Heptane_Prescribed_2cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR NIST_Heptane_Prescribed_1cm.fds
$QFDS $DEBUG -p 32  $QUEUE -d $INDIR NIST_Heptane_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR NIST_Methanol_Prescribed_2cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR NIST_Methanol_Prescribed_1cm.fds
$QFDS $DEBUG -p 32  $QUEUE -d $INDIR NIST_Methanol_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Methane_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Methane_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Methane_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_20kW_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_20kW_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Propane_20kW_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_34kW_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_34kW_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Propane_34kW_Prescribed_0p5cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_50kW_Prescribed_2cm.fds
$QFDS $DEBUG -p 20  $QUEUE -d $INDIR NIST_Propane_50kW_Prescribed_1cm.fds
$QFDS $DEBUG -p 40  $QUEUE -d $INDIR NIST_Propane_50kW_Prescribed_0p5cm.fds
 
$QFDS $DEBUG -p 108 $QUEUE -d $INDIR NIST_Methanol_1m_pan_1cm_grid.fds
$QFDS $DEBUG -p 24  $QUEUE -d $INDIR NIST_Methanol_1m_pan_2cm_grid.fds
$QFDS $DEBUG -p 10  $QUEUE -d $INDIR NIST_Methanol_1m_pan_4cm_grid.fds
$QFDS $DEBUG -p 108 $QUEUE -d $INDIR NIST_Methanol_1m_pan_1cm_grid_predicted.fds
$QFDS $DEBUG -p 24  $QUEUE -d $INDIR NIST_Methanol_1m_pan_2cm_grid_predicted.fds
$QFDS $DEBUG -p 10  $QUEUE -d $INDIR NIST_Methanol_1m_pan_4cm_grid_predicted.fds

echo FDS cases submitted
