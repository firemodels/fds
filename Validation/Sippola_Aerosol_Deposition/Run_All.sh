#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Sippola_Test_01
$RUNFDS $INDIR Sippola_Test_02
$RUNFDS $INDIR Sippola_Test_03
$RUNFDS $INDIR Sippola_Test_04
$RUNFDS $INDIR Sippola_Test_05
$RUNFDS $INDIR Sippola_Test_06
$RUNFDS $INDIR Sippola_Test_07
$RUNFDS $INDIR Sippola_Test_08
$RUNFDS $INDIR Sippola_Test_09
$RUNFDS $INDIR Sippola_Test_10
$RUNFDS $INDIR Sippola_Test_11
$RUNFDS $INDIR Sippola_Test_12
$RUNFDS $INDIR Sippola_Test_13
$RUNFDS $INDIR Sippola_Test_14
$RUNFDS $INDIR Sippola_Test_15
$RUNFDS $INDIR Sippola_Test_16

echo FDS cases submitted
