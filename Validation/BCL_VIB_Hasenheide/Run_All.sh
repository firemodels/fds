#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S1_R1.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S1_R2.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S1_R3.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S2_R1.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S2_R2.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S2_R3.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S3_R1.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S3_R2.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S4_R1.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_160_S4_R2.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_220_S1_R1.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_220_S1_R2.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_220_S1_R3.fds
$QFDS $DEBUG -p 24 $QUEUE -d $INDIR  Hasenheide_Hep_220_S1_R4.fds
