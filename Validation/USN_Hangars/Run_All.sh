#!/bin/bash

# This script runs the a set of Validation Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/../..
export BASEDIR=`pwd`
export INDIR=Current_Results
source ../../FDS_Compilation/SET_MYFDSENV.sh intel64 run

# to override FDSMPI defined in above script, remove comment
# from  line below and define your own FDSMPI location
#export FDSMPI=override FDSMPI defined in above script

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

