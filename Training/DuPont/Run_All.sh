#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=~/FDS/FDS5/bin/fds5_linux_64
#export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR ASTM_E_648_12 
$RUNFDS $INDIR ASTM_E_648_13
$RUNFDS $INDIR ASTM_E_648_14 
$RUNFDS $INDIR ASTM_E_648_15
$RUNFDS $INDIR ASTM_E_648_16
$RUNFDS $INDIR ASTM_E_648_17 

echo FDS cases submitted
