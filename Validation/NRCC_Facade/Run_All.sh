#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR NRCC_Facade_Win_1_05_MW 
$RUNFDS $INDIR NRCC_Facade_Win_1_06_MW 
$RUNFDS $INDIR NRCC_Facade_Win_1_08_MW 
$RUNFDS $INDIR NRCC_Facade_Win_2_05_MW 
$RUNFDS $INDIR NRCC_Facade_Win_2_06_MW 
$RUNFDS $INDIR NRCC_Facade_Win_2_08_MW 
$RUNFDS $INDIR NRCC_Facade_Win_2_10_MW 
$RUNFDS $INDIR NRCC_Facade_Win_3_05_MW 
$RUNFDS $INDIR NRCC_Facade_Win_3_06_MW 
$RUNFDS $INDIR NRCC_Facade_Win_3_08_MW 
$RUNFDS $INDIR NRCC_Facade_Win_3_10_MW 
$RUNFDS $INDIR NRCC_Facade_Win_4_05_MW 
$RUNFDS $INDIR NRCC_Facade_Win_4_06_MW 
$RUNFDS $INDIR NRCC_Facade_Win_4_08_MW 
$RUNFDS $INDIR NRCC_Facade_Win_4_10_MW 
$RUNFDS $INDIR NRCC_Facade_Win_5_05_MW 
$RUNFDS $INDIR NRCC_Facade_Win_5_06_MW 
$RUNFDS $INDIR NRCC_Facade_Win_5_08_MW 
$RUNFDS $INDIR NRCC_Facade_Win_5_10_MW 

echo FDS cases submitted
