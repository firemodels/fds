#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

# good burn: 136, 114, marginal burn: 41,80 no burn: 30,60

#$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn136_120D_30S_35L.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn114_120D_18S_35L.fds
#$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn41_60D_24S_30L.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn80_120D_3S_40L.fds
#$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn30_30D_3S_30L.fds
$QFDS $DEBUG $QUEUE -p 16 -d $INDIR burn60_120D_0S_30L.fds

echo FDS cases submitted
