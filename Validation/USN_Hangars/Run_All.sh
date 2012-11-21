#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export FDSMPI=$SVNROOT/FDS_Compilation/mpi_intel_linux_64/fds_mpi_intel_linux_64
export RUNFDSMPI=$SVNROOT/Utilities/Scripts/runfdsmpi.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_01
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_02
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_03
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_04
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_05
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_06
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_07
$RUNFDSMPI 5 $INDIR USN_Hawaii_Test_11

$RUNFDSMPI 5 $INDIR USN_Iceland_Test_01
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_02
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_03
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_04
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_05
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_06
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_07
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_09
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_10
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_11
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_12
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_13
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_14
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_15
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_17
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_18
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_19
$RUNFDSMPI 5 $INDIR USN_Iceland_Test_20

