#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_tga_N2_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_tga_N2_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_tga_air_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_tga_air_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_mcc_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_S701_mcc_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_I701_tga_N2_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_I701_tga_N2_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_I701_mcc_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_I701_mcc_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_25_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_25_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_50_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_50_v2.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_75_v1.fds
$QFDS $DEBUG $QUEUE -d $INDIR CHRISTIFIRE_C701_cone_75_v2.fds

echo FDS cases submitted
