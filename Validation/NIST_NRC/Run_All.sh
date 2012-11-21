#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR NIST_NRC_01 
$RUNFDS $INDIR NIST_NRC_02 
$RUNFDS $INDIR NIST_NRC_03 
$RUNFDS $INDIR NIST_NRC_04 
$RUNFDS $INDIR NIST_NRC_05 
$RUNFDS $INDIR NIST_NRC_07 
$RUNFDS $INDIR NIST_NRC_08 
$RUNFDS $INDIR NIST_NRC_09 
$RUNFDS $INDIR NIST_NRC_10 
$RUNFDS $INDIR NIST_NRC_13 
$RUNFDS $INDIR NIST_NRC_14 
$RUNFDS $INDIR NIST_NRC_15 
$RUNFDS $INDIR NIST_NRC_16 
$RUNFDS $INDIR NIST_NRC_17 
$RUNFDS $INDIR NIST_NRC_18 

echo FDS cases submitted
