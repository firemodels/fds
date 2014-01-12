#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_cold_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_hot_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum _crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_noplenum_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 200_once_plenum_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_cold_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_hot_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum _crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_noplenum_server9_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle1.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle2.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle3.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle4.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle5.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_aisle6.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_crah.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_server7_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_server7_cable.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_server9_board.fds
$QFDS $DEBUG -p 4 -r $QUEUE -d $INDIR 50_once_plenum_server9_cable.fds
 
echo FDS cases submitted
