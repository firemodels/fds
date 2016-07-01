#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISONG3.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept4.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept5.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept8.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept9.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISONylon10.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOPP11.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHeptD12.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHeptD13.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOPropD14.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOProp15.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOStyrene16.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOStyrene17.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOPP18.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept19.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOToluene20.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOStyrene21.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept22.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept23.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept24.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept25.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept26.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept27.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOHept28.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOToluene29.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISOPropanol30.fds
$QFDS $DEBUG $QUEUE -p 2 -d $INDIR ISONG32.fds

echo FDS cases submitted
