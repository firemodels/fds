#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

#$QFDS -p 1 -n 12 $DEBUG $QUEUE -d $INDIR NFRL_8MW_Burner.fds
#$QFDS -p 1 -n 12 $DEBUG $QUEUE -d $INDIR NFRL_20MW_Burner.fds

echo FDS cases submitted
