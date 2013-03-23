#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR UL_NIST_Vents_Test_1
$RUNFDS $INDIR UL_NIST_Vents_Test_2
$RUNFDS $INDIR UL_NIST_Vents_Test_3
$RUNFDS $INDIR UL_NIST_Vents_Test_4

echo FDS cases submitted
