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

$QFDS -r $qq -d $INDIR moody_dpdx=-0p01_N16.fds  
$QFDS -r $qq -d $INDIR moody_dpdx=-0p01_N32.fds  
$QFDS -r $qq -d $INDIR moody_dpdx=-0p01_N8.fds   
$QFDS -r $qq -d $INDIR moody_dpdx=-100_N16.fds   
$QFDS -r $qq -d $INDIR moody_dpdx=-100_N32.fds   
$QFDS -r $qq -d $INDIR moody_dpdx=-100_N8.fds    
$QFDS -r $qq -d $INDIR moody_dpdx=-1_N16.fds     
$QFDS -r $qq -d $INDIR moody_dpdx=-1_N32.fds     
$QFDS -r $qq -d $INDIR moody_dpdx=-1_N8.fds      
$QFDS -r $qq -d $INDIR poiseuille_N16_mu025.fds  
$QFDS -r $qq -d $INDIR poiseuille_N32_mu025.fds  
$QFDS -r $qq -d $INDIR poiseuille_N64_mu0125.fds 
$QFDS -r $qq -d $INDIR poiseuille_N64_mu025.fds  
$QFDS -r $qq -d $INDIR poiseuille_N8_mu025.fds   

$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-1_N8.fds     
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-1_N8.fds      
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-1_N8.fds       
$QFDS -r $qq -d $INDIR z0=p1_dpdx=-1_N8.fds        
$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-1_N16.fds    
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-1_N16.fds     
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-1_N16.fds      
$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-100_N8.fds   
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-100_N8.fds    
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-100_N8.fds     
$QFDS -r $qq -d $INDIR z0=p1_dpdx=-100_N8.fds      
$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-100_N16.fds  
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-100_N16.fds   
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-100_N16.fds    
$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-p0001_N8.fds 
$QFDS -r $qq -d $INDIR z0=p0001_dpdx=-p01_N8.fds   
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-p0001_N8.fds  
$QFDS -r $qq -d $INDIR z0=p001_dpdx=-p01_N8.fds    
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-p0001_N8.fds   
$QFDS -r $qq -d $INDIR z0=p01_dpdx=-p01_N8.fds     
$QFDS -r $qq -d $INDIR z0=p1_dpdx=-p0001_N8.fds    
$QFDS -r $qq -d $INDIR z0=p1_dpdx=-p01_N8.fds      
$QFDS -r $qq -d $INDIR z0=p02_dpdx=-1_N16.fds      
$QFDS -r $qq -d $INDIR z0=p02_dpdx=-100_N16.fds    

echo FDS cases submitted
