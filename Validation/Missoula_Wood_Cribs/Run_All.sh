#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_1.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_2.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_3.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_4.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_5.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_6.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_7.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_8a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_8b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_8c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_8d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_8e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_9.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_10.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_11.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_12.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_13.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_14.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_15a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_15b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_15c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_15d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_15e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_16a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_16c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_16d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_16e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_17.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_18a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_18b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_18c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_18d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_18e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_19a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_19c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_19d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_19e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_20a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_20b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_20c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_20d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_20e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_21a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_21b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_21c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_21d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_21e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_22a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_22c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_22d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_22e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_23.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_24a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_24c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_24d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_24e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_25.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_26.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_27.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_28a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_28c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_28d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_28e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_29a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_29b.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_29c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_29d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_29e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_30.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_31a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_31c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_31d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_31e.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_32.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_33a.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_33c.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_33d.fds
$QFDS $DEBUG $QUEUE -p 1 -d $INDIR missoula_cribs_33e.fds

echo FDS cases submitted

