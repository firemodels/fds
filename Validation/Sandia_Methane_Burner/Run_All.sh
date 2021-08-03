#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner01.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner02.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner03.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner04.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner05.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner06.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner07.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner08.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner09.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner10.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner11.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner12.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner13.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner14.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner15.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner16.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner17.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner18.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner19.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner20.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner21.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner22.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner23.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner24.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner25.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner26.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner27.fds
$QFDS $DEBUG $QUEUE -p 12 -d $INDIR Burner28.fds

echo FDS cases submitted
