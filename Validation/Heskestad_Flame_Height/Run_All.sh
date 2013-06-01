#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR="Current_Results"
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR Qs=10000_RI=05.fds
$QFDS -r $qq -d $INDIR Qs=1000_RI=05.fds  
$QFDS -r $qq -d $INDIR Qs=100_RI=05.fds   
$QFDS -r $qq -d $INDIR Qs=10_RI=05.fds    
$QFDS -r $qq -d $INDIR Qs=1_RI=05.fds     
$QFDS -r $qq -d $INDIR Qs=2000_RI=05.fds  
$QFDS -r $qq -d $INDIR Qs=200_RI=05.fds   
$QFDS -r $qq -d $INDIR Qs=20_RI=05.fds    
$QFDS -r $qq -d $INDIR Qs=2_RI=05.fds     
$QFDS -r $qq -d $INDIR Qs=5000_RI=05.fds  
$QFDS -r $qq -d $INDIR Qs=500_RI=05.fds   
$QFDS -r $qq -d $INDIR Qs=50_RI=05.fds    
$QFDS -r $qq -d $INDIR Qs=5_RI=05.fds     
$QFDS -r $qq -d $INDIR Qs=p1_RI=05.fds    
$QFDS -r $qq -d $INDIR Qs=p2_RI=05.fds    
$QFDS -r $qq -d $INDIR Qs=p5_RI=05.fds    

$QFDS -r $qq -d $INDIR Qs=10000_RI=10.fds 
$QFDS -r $qq -d $INDIR Qs=1000_RI=10.fds  
$QFDS -r $qq -d $INDIR Qs=100_RI=10.fds   
$QFDS -r $qq -d $INDIR Qs=10_RI=10.fds    
$QFDS -r $qq -d $INDIR Qs=1_RI=10.fds     
$QFDS -r $qq -d $INDIR Qs=2000_RI=10.fds  
$QFDS -r $qq -d $INDIR Qs=200_RI=10.fds   
$QFDS -r $qq -d $INDIR Qs=20_RI=10.fds    
$QFDS -r $qq -d $INDIR Qs=2_RI=10.fds     
$QFDS -r $qq -d $INDIR Qs=5000_RI=10.fds  
$QFDS -r $qq -d $INDIR Qs=500_RI=10.fds   
$QFDS -r $qq -d $INDIR Qs=50_RI=10.fds    
$QFDS -r $qq -d $INDIR Qs=5_RI=10.fds     
$QFDS -r $qq -d $INDIR Qs=p1_RI=10.fds    
$QFDS -r $qq -d $INDIR Qs=p2_RI=10.fds    
$QFDS -r $qq -d $INDIR Qs=p5_RI=10.fds    

$QFDS -r $qq -d $INDIR Qs=10000_RI=20.fds 
$QFDS -r $qq -d $INDIR Qs=1000_RI=20.fds  
$QFDS -r $qq -d $INDIR Qs=100_RI=20.fds   
$QFDS -r $qq -d $INDIR Qs=10_RI=20.fds    
$QFDS -r $qq -d $INDIR Qs=1_RI=20.fds     
$QFDS -r $qq -d $INDIR Qs=2000_RI=20.fds  
$QFDS -r $qq -d $INDIR Qs=200_RI=20.fds   
$QFDS -r $qq -d $INDIR Qs=20_RI=20.fds    
$QFDS -r $qq -d $INDIR Qs=2_RI=20.fds     
$QFDS -r $qq -d $INDIR Qs=5000_RI=20.fds  
$QFDS -r $qq -d $INDIR Qs=500_RI=20.fds   
$QFDS -r $qq -d $INDIR Qs=50_RI=20.fds    
$QFDS -r $qq -d $INDIR Qs=5_RI=20.fds     
$QFDS -r $qq -d $INDIR Qs=p1_RI=20.fds    
$QFDS -r $qq -d $INDIR Qs=p2_RI=20.fds    
$QFDS -r $qq -d $INDIR Qs=p5_RI=20.fds    

echo FDS cases submitted
