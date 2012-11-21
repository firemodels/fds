#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Vettori_OBSTRUCTED_CORNER_FAST 
$RUNFDS $INDIR Vettori_OBSTRUCTED_CORNER_MED  
$RUNFDS $INDIR Vettori_OBSTRUCTED_CORNER_SLOW 
$RUNFDS $INDIR Vettori_OBSTRUCTED_OPEN_FAST   
$RUNFDS $INDIR Vettori_OBSTRUCTED_OPEN_MED    
$RUNFDS $INDIR Vettori_OBSTRUCTED_OPEN_SLOW   
$RUNFDS $INDIR Vettori_OBSTRUCTED_WALL_FAST   
$RUNFDS $INDIR Vettori_OBSTRUCTED_WALL_MED    
$RUNFDS $INDIR Vettori_OBSTRUCTED_WALL_SLOW   
$RUNFDS $INDIR Vettori_SMOOTH_CORNER_FAST     
$RUNFDS $INDIR Vettori_SMOOTH_CORNER_MED      
$RUNFDS $INDIR Vettori_SMOOTH_CORNER_SLOW     
$RUNFDS $INDIR Vettori_SMOOTH_OPEN_FAST       
$RUNFDS $INDIR Vettori_SMOOTH_OPEN_MED        
$RUNFDS $INDIR Vettori_SMOOTH_OPEN_SLOW       
$RUNFDS $INDIR Vettori_SMOOTH_WALL_FAST       
$RUNFDS $INDIR Vettori_SMOOTH_WALL_MED        
$RUNFDS $INDIR Vettori_SMOOTH_WALL_SLOW       

echo FDS cases submitted
