#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_1_Lower.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_1_Upper.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_2_Lower.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_2_Upper.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_3.fds      
$QFDS $DEBUG $QUEUE -d $INDIR PRISME_LK_4.fds      

$QFDS $DEBUG $QUEUE -d $INDIR PRS_D1.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_D2.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_D3.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_D4.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_D5.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_D6.fds

$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D1.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D2.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D3.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D4.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D5a.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D5.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D6a.fds
$QFDS $DEBUG $QUEUE -d $INDIR PRS_SI_D6.fds

$QFDS $DEBUG $QUEUE -p 16 -d $INDIR PR2_FES_PA2a.fds

echo FDS cases submitted
