#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_CORNER_FAST.fds 
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_CORNER_MED.fds  
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_CORNER_SLOW.fds 
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_OPEN_FAST.fds   
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_OPEN_MED.fds    
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_OPEN_SLOW.fds   
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_WALL_FAST.fds   
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_WALL_MED.fds    
$QFDS -r $qq -d $INDIR Vettori_OBSTRUCTED_WALL_SLOW.fds   
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_CORNER_FAST.fds     
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_CORNER_MED.fds      
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_CORNER_SLOW.fds     
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_OPEN_FAST.fds       
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_OPEN_MED.fds        
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_OPEN_SLOW.fds       
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_WALL_FAST.fds       
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_WALL_MED.fds        
$QFDS -r $qq -d $INDIR Vettori_SMOOTH_WALL_SLOW.fds      

echo FDS cases submitted
