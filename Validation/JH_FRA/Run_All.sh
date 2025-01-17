#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_01.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR JH_FRA_compartment_02.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_03.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_03A.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_04.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_11.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR JH_FRA_compartment_12.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_13.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_14.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_21.fds
$QFDS $DEBUG $QUEUE -p 4 -d $INDIR JH_FRA_compartment_22.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_23.fds
$QFDS $DEBUG $QUEUE -p 6 -d $INDIR JH_FRA_compartment_24.fds

echo FDS cases submitted
