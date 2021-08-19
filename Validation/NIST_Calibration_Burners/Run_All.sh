#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

#$QFDS -p 424 $DEBUG $QUEUE -d $INDIR NIST_8MW_Burner.fds
#$QFDS -p 400 $DEBUG $QUEUE -d $INDIR NIST_20MW_Burner.fds
#$QFDS -p 240 $DEBUG $QUEUE -d $INDIR NIST_20MW_Ports.fds

echo FDS cases submitted
