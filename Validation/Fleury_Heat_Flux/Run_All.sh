#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Fleury_1t1_100_kW
$RUNFDS $INDIR Fleury_1t1_150_kW
$RUNFDS $INDIR Fleury_1t1_200_kW
$RUNFDS $INDIR Fleury_1t1_250_kW
$RUNFDS $INDIR Fleury_1t1_300_kW
$RUNFDS $INDIR Fleury_2t1_100_kW
$RUNFDS $INDIR Fleury_2t1_150_kW
$RUNFDS $INDIR Fleury_2t1_200_kW
$RUNFDS $INDIR Fleury_2t1_250_kW
$RUNFDS $INDIR Fleury_2t1_300_kW
$RUNFDS $INDIR Fleury_3t1_100_kW
$RUNFDS $INDIR Fleury_3t1_150_kW
$RUNFDS $INDIR Fleury_3t1_200_kW
$RUNFDS $INDIR Fleury_3t1_250_kW
$RUNFDS $INDIR Fleury_3t1_300_kW

echo FDS cases submitted
