#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_CORNER_FAST.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_CORNER_MED.fds  
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_CORNER_SLOW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_OPEN_FAST.fds   
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_OPEN_MED.fds    
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_OPEN_SLOW.fds   
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_WALL_FAST.fds   
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_WALL_MED.fds    
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_OBSTRUCTED_WALL_SLOW.fds   
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_CORNER_FAST.fds     
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_CORNER_MED.fds      
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_CORNER_SLOW.fds     
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_OPEN_FAST.fds       
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_OPEN_MED.fds        
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_OPEN_SLOW.fds       
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_WALL_FAST.fds       
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_WALL_MED.fds        
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_SMOOTH_WALL_SLOW.fds      

echo FDS cases submitted
