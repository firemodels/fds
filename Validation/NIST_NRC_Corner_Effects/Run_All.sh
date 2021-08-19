#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 15 -d $INDIR corner_200_kW.fds
$QFDS $DEBUG $QUEUE -p 15 -d $INDIR corner_300_kW.fds
$QFDS $DEBUG $QUEUE -p 15 -d $INDIR corner_400_kW.fds
$QFDS $DEBUG $QUEUE -p 17 -d $INDIR wall_200_kW.fds
$QFDS $DEBUG $QUEUE -p 17 -d $INDIR wall_300_kW.fds
$QFDS $DEBUG $QUEUE -p 17 -d $INDIR wall_400_kW.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_01.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_02.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_03.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_04.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_05.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_06.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_07.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_08.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_09.fds
$QFDS $DEBUG $QUEUE -p 6  -d $INDIR cabinet_10.fds
$QFDS $DEBUG $QUEUE -p 13 -d $INDIR cabinet_11.fds
$QFDS $DEBUG $QUEUE -p 13 -d $INDIR cabinet_12.fds

