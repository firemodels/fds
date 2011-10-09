#!/bin/bash -f

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

$RUNFDSMPI 5 $INDIR FM_SNL_01
$RUNFDSMPI 5 $INDIR FM_SNL_02
$RUNFDSMPI 5 $INDIR FM_SNL_03
$RUNFDSMPI 5 $INDIR FM_SNL_04
$RUNFDSMPI 5 $INDIR FM_SNL_05
$RUNFDSMPI 5 $INDIR FM_SNL_06
$RUNFDSMPI 5 $INDIR FM_SNL_07
$RUNFDSMPI 5 $INDIR FM_SNL_08
$RUNFDSMPI 5 $INDIR FM_SNL_09
$RUNFDSMPI 5 $INDIR FM_SNL_10
$RUNFDSMPI 5 $INDIR FM_SNL_11
$RUNFDSMPI 5 $INDIR FM_SNL_12
$RUNFDSMPI 5 $INDIR FM_SNL_13
$RUNFDSMPI 5 $INDIR FM_SNL_14
$RUNFDSMPI 5 $INDIR FM_SNL_15
$RUNFDSMPI 5 $INDIR FM_SNL_16
$RUNFDSMPI 5 $INDIR FM_SNL_17
$RUNFDSMPI 5 $INDIR FM_SNL_21
$RUNFDSMPI 5 $INDIR FM_SNL_22
