#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh -q serial -y

export RESOURCE_MANAGER='SLURM'

$QFDS $DEBUG -w 48:00:00 -q serial  -p 3 -n 3 -d $INDIR TUS_facade_I-1-2.fds
$QFDS $DEBUG -w 48:00:00 -q serial  -p 3 -n 3 -d $INDIR TUS_facade_1-1-2-3.fds
$QFDS $DEBUG -w 48:00:00 -q serial  -p 3 -n 3 -d $INDIR TUS_facade_2-1-2-3.fds
