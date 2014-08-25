#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR Qs=10000_RI=05.fds
$QFDS  $QUEUE -d $INDIR Qs=1000_RI=05.fds  
$QFDS  $QUEUE -d $INDIR Qs=100_RI=05.fds   
$QFDS  $QUEUE -d $INDIR Qs=10_RI=05.fds    
$QFDS  $QUEUE -d $INDIR Qs=1_RI=05.fds     
$QFDS  $QUEUE -d $INDIR Qs=2000_RI=05.fds  
$QFDS  $QUEUE -d $INDIR Qs=200_RI=05.fds   
$QFDS  $QUEUE -d $INDIR Qs=20_RI=05.fds    
$QFDS  $QUEUE -d $INDIR Qs=2_RI=05.fds     
$QFDS  $QUEUE -d $INDIR Qs=5000_RI=05.fds  
$QFDS  $QUEUE -d $INDIR Qs=500_RI=05.fds   
$QFDS  $QUEUE -d $INDIR Qs=50_RI=05.fds    
$QFDS  $QUEUE -d $INDIR Qs=5_RI=05.fds     
$QFDS  $QUEUE -d $INDIR Qs=p1_RI=05.fds    
$QFDS  $QUEUE -d $INDIR Qs=p2_RI=05.fds    
$QFDS  $QUEUE -d $INDIR Qs=p5_RI=05.fds    

$QFDS  $QUEUE -d $INDIR Qs=10000_RI=10.fds 
$QFDS  $QUEUE -d $INDIR Qs=1000_RI=10.fds  
$QFDS  $QUEUE -d $INDIR Qs=100_RI=10.fds   
$QFDS  $QUEUE -d $INDIR Qs=10_RI=10.fds    
$QFDS  $QUEUE -d $INDIR Qs=1_RI=10.fds     
$QFDS  $QUEUE -d $INDIR Qs=2000_RI=10.fds  
$QFDS  $QUEUE -d $INDIR Qs=200_RI=10.fds   
$QFDS  $QUEUE -d $INDIR Qs=20_RI=10.fds    
$QFDS  $QUEUE -d $INDIR Qs=2_RI=10.fds     
$QFDS  $QUEUE -d $INDIR Qs=5000_RI=10.fds  
$QFDS  $QUEUE -d $INDIR Qs=500_RI=10.fds   
$QFDS  $QUEUE -d $INDIR Qs=50_RI=10.fds    
$QFDS  $QUEUE -d $INDIR Qs=5_RI=10.fds     
$QFDS  $QUEUE -d $INDIR Qs=p1_RI=10.fds    
$QFDS  $QUEUE -d $INDIR Qs=p2_RI=10.fds    
$QFDS  $QUEUE -d $INDIR Qs=p5_RI=10.fds    

$QFDS  $QUEUE -d $INDIR Qs=10000_RI=20.fds 
$QFDS  $QUEUE -d $INDIR Qs=1000_RI=20.fds  
$QFDS  $QUEUE -d $INDIR Qs=100_RI=20.fds   
$QFDS  $QUEUE -d $INDIR Qs=10_RI=20.fds    
$QFDS  $QUEUE -d $INDIR Qs=1_RI=20.fds     
$QFDS  $QUEUE -d $INDIR Qs=2000_RI=20.fds  
$QFDS  $QUEUE -d $INDIR Qs=200_RI=20.fds   
$QFDS  $QUEUE -d $INDIR Qs=20_RI=20.fds    
$QFDS  $QUEUE -d $INDIR Qs=2_RI=20.fds     
$QFDS  $QUEUE -d $INDIR Qs=5000_RI=20.fds  
$QFDS  $QUEUE -d $INDIR Qs=500_RI=20.fds   
$QFDS  $QUEUE -d $INDIR Qs=50_RI=20.fds    
$QFDS  $QUEUE -d $INDIR Qs=5_RI=20.fds     
$QFDS  $QUEUE -d $INDIR Qs=p1_RI=20.fds    
$QFDS  $QUEUE -d $INDIR Qs=p2_RI=20.fds    
$QFDS  $QUEUE -d $INDIR Qs=p5_RI=20.fds    

echo FDS cases submitted
