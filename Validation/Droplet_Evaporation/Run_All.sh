#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR fujita_u2p0_hp0.fds
$QFDS $DEBUG $QUEUE -d $INDIR fujita_u2p0_hp3.fds
$QFDS $DEBUG $QUEUE -d $INDIR fujita_up8_hp0.fds
$QFDS $DEBUG $QUEUE -d $INDIR fujita_up8_hp3.fds
$QFDS $DEBUG $QUEUE -d $INDIR gavin_557.fds
$QFDS $DEBUG $QUEUE -d $INDIR gavin_769.fds
$QFDS $DEBUG $QUEUE -d $INDIR kolaitis_decane.fds
$QFDS $DEBUG $QUEUE -d $INDIR kolaitis_ethanol.fds
$QFDS $DEBUG $QUEUE -d $INDIR kolaitis_heptane_1.fds
$QFDS $DEBUG $QUEUE -d $INDIR kolaitis_heptane_2.fds
$QFDS $DEBUG $QUEUE -d $INDIR maqua_acetone.fds
$QFDS $DEBUG $QUEUE -d $INDIR maqua_ethanol.fds
$QFDS $DEBUG $QUEUE -d $INDIR taflin_43p9.fds
$QFDS $DEBUG $QUEUE -d $INDIR taflin_56p6.fds