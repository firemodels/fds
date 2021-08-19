#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-100-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-100-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-100-p39c.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-100-p39d.fds

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-100-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-100-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-100-p39c.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-100-p39d.fds

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-100-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-100-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-100-p39c.fds

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-200-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-200-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-200-p39c.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-2p5-200-p39d.fds

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-200-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-200-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-200-p39c.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-5p0-200-p39d.fds

$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-200-p39a.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-200-p39b.fds
$QFDS -p 6  $DEBUG $QUEUE -d $INDIR NIST_SDG-10p0-200-p39c.fds

echo FDS cases submitted







