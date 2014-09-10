#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR 7p1_cm_methane_4mm.fds
#$QFDS $DEBUG $QUEUE -d $INDIR 7p1_cm_methane_2mm.fds
$QFDS $DEBUG -p 16 $QUEUE -d $INDIR 7p1_cm_methane_2mm_16mesh.fds
#$QFDS $DEBUG -p 128 $QUEUE -d $INDIR 7p1_cm_methane_1mm_128mesh_les.fds
#$QFDS $DEBUG -p 128 $QUEUE -d $INDIR 7p1_cm_methane_1mm_128mesh_dns.fds

echo FDS cases submitted
