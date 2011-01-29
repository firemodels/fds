#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR NIST_RSE_1994_100 
$RUNFDS $INDIR NIST_RSE_1994_150 
$RUNFDS $INDIR NIST_RSE_1994_200 
$RUNFDS $INDIR NIST_RSE_1994_300 
$RUNFDS $INDIR NIST_RSE_1994_400 
$RUNFDS $INDIR NIST_RSE_1994_500 
$RUNFDS $INDIR NIST_RSE_1994_50  
$RUNFDS $INDIR NIST_RSE_1994_600 
$RUNFDS $INDIR NIST_RSE_1994_75  

echo FDS cases submitted
