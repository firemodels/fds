#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 1   $QUEUE -d $INDIR Sandia_He_1m_dx20cm.fds
$QFDS $DEBUG -p 1   $QUEUE -d $INDIR Sandia_He_1m_dx10cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Sandia_He_1m_dx6cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Sandia_He_1m_dx3cm.fds
$QFDS $DEBUG -p 64  $QUEUE -d $INDIR Sandia_He_1m_dx1p5cm.fds

$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test14_dx6cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test17_dx6cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test24_dx6cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_H2_1m_Test35_dx6cm.fds

$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test14_dx3cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test17_dx3cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test24_dx3cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_H2_1m_Test35_dx3cm.fds

$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test14_dx1p5cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test17_dx1p5cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_CH4_1m_Test24_dx1p5cm.fds
$QFDS $DEBUG -p 16  $QUEUE -d $INDIR Sandia_H2_1m_Test35_dx1p5cm.fds

echo FDS cases submitted
