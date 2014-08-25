#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR ISONG3.fds
$QFDS  $QUEUE -d $INDIR ISOHept4.fds
$QFDS  $QUEUE -d $INDIR ISOHept5.fds
$QFDS  $QUEUE -d $INDIR ISOHept8.fds
$QFDS  $QUEUE -d $INDIR ISOHept9.fds
$QFDS  $QUEUE -d $INDIR ISONylon10.fds
$QFDS  $QUEUE -d $INDIR ISOPP11.fds
$QFDS  $QUEUE -d $INDIR ISOHeptD12.fds
$QFDS  $QUEUE -d $INDIR ISOHeptD13.fds
$QFDS  $QUEUE -d $INDIR ISOPropD14.fds
$QFDS  $QUEUE -d $INDIR ISOProp15.fds
$QFDS  $QUEUE -d $INDIR ISOStyrene16.fds
$QFDS  $QUEUE -d $INDIR ISOStyrene17.fds
$QFDS  $QUEUE -d $INDIR ISOPP18.fds
$QFDS  $QUEUE -d $INDIR ISOHept19.fds
$QFDS  $QUEUE -d $INDIR ISOToluene20.fds
$QFDS  $QUEUE -d $INDIR ISOStyrene21.fds
$QFDS  $QUEUE -d $INDIR ISOHept22.fds
$QFDS  $QUEUE -d $INDIR ISOHept23.fds
$QFDS  $QUEUE -d $INDIR ISOHept24.fds
$QFDS  $QUEUE -d $INDIR ISOHept25.fds
$QFDS  $QUEUE -d $INDIR ISOHept26.fds
$QFDS  $QUEUE -d $INDIR ISOHept27.fds
$QFDS  $QUEUE -d $INDIR ISOHept28.fds
$QFDS  $QUEUE -d $INDIR ISOToluene29.fds
$QFDS  $QUEUE -d $INDIR ISOPropanol30.fds
$QFDS  $QUEUE -d $INDIR ISONG32.fds

echo FDS cases submitted
