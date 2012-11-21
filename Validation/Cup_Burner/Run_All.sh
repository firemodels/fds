#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Cup_C7H16_Ar
$RUNFDS $INDIR Cup_C7H16_CO2
$RUNFDS $INDIR Cup_C7H16_He
$RUNFDS $INDIR Cup_C7H16_N2
$RUNFDS $INDIR Cup_CH4_Ar
$RUNFDS $INDIR Cup_CH4_CO2
$RUNFDS $INDIR Cup_CH4_He
$RUNFDS $INDIR Cup_CH4_N2

echo FDS cases submitted
