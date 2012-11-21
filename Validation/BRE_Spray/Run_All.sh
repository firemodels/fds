#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR BRE_Spray_A_1 
$RUNFDS $INDIR BRE_Spray_A_2
$RUNFDS $INDIR BRE_Spray_A_3 
$RUNFDS $INDIR BRE_Spray_A_4 
$RUNFDS $INDIR BRE_Spray_A_5 
$RUNFDS $INDIR BRE_Spray_A_6 
$RUNFDS $INDIR BRE_Spray_A_7 
$RUNFDS $INDIR BRE_Spray_A_8
$RUNFDS $INDIR BRE_Spray_B_1 
$RUNFDS $INDIR BRE_Spray_B_2
$RUNFDS $INDIR BRE_Spray_B_3 
$RUNFDS $INDIR BRE_Spray_B_4 
$RUNFDS $INDIR BRE_Spray_B_5 
$RUNFDS $INDIR BRE_Spray_B_6 
$RUNFDS $INDIR BRE_Spray_B_7 
$RUNFDS $INDIR BRE_Spray_B_8
$RUNFDS $INDIR BRE_Spray_D_1 
$RUNFDS $INDIR BRE_Spray_D_2
$RUNFDS $INDIR BRE_Spray_D_3 
$RUNFDS $INDIR BRE_Spray_D_4 
$RUNFDS $INDIR BRE_Spray_D_5 
$RUNFDS $INDIR BRE_Spray_D_6 
$RUNFDS $INDIR BRE_Spray_D_7 
$RUNFDS $INDIR BRE_Spray_D_8
 
echo FDS cases submitted
