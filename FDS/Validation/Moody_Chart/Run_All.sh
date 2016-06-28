#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N16.fds  
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N32.fds  
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-0p01_N8.fds   
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N16.fds   
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N32.fds   
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-100_N8.fds    
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N16.fds     
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N32.fds     
$QFDS $DEBUG $QUEUE -d $INDIR moody_dpdx=-1_N8.fds      
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N16_mu025.fds  
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N32_mu025.fds  
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N64_mu0125.fds 
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N64_mu025.fds  
$QFDS $DEBUG $QUEUE -d $INDIR poiseuille_N8_mu025.fds   

$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-1_N8.fds     
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-1_N8.fds      
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-1_N8.fds       
$QFDS $DEBUG $QUEUE -d $INDIR z0=p1_dpdx=-1_N8.fds        
$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-1_N16.fds    
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-1_N16.fds     
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-1_N16.fds      
$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-100_N8.fds   
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-100_N8.fds    
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-100_N8.fds     
$QFDS $DEBUG $QUEUE -d $INDIR z0=p1_dpdx=-100_N8.fds      
$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-100_N16.fds  
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-100_N16.fds   
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-100_N16.fds    
$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-p0001_N8.fds 
$QFDS $DEBUG $QUEUE -d $INDIR z0=p0001_dpdx=-p01_N8.fds   
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-p0001_N8.fds  
$QFDS $DEBUG $QUEUE -d $INDIR z0=p001_dpdx=-p01_N8.fds    
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-p0001_N8.fds   
$QFDS $DEBUG $QUEUE -d $INDIR z0=p01_dpdx=-p01_N8.fds     
$QFDS $DEBUG $QUEUE -d $INDIR z0=p1_dpdx=-p0001_N8.fds    
$QFDS $DEBUG $QUEUE -d $INDIR z0=p1_dpdx=-p01_N8.fds      
$QFDS $DEBUG $QUEUE -d $INDIR z0=p02_dpdx=-1_N16.fds      
$QFDS $DEBUG $QUEUE -d $INDIR z0=p02_dpdx=-100_N16.fds    

echo FDS cases submitted
