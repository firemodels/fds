#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Steckler_010.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_011.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_012.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_013.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_014.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_016.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_017.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_018.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_019.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_020.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_021.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_022.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_023.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_030.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_041.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_114.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_116.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_122.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_144.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_160.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_161.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_162.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_163.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_164.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_165.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_166.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_167.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_210.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_212.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_220.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_221.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_224.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_240.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_242.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_310.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_324.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_410.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_510.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_512.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_513.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_514.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_517.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_520.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_521.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_522.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_524.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_540.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_541.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_542.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_544.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_610.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_612.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_622.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_710.fds
$QFDS $DEBUG $QUEUE -d $INDIR Steckler_810.fds

echo FDS cases submitted
