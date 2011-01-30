#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR moody_dpdx=-0p01_N16  
$RUNFDS $INDIR moody_dpdx=-0p01_N32  
$RUNFDS $INDIR moody_dpdx=-0p01_N8   
$RUNFDS $INDIR moody_dpdx=-100_N16   
$RUNFDS $INDIR moody_dpdx=-100_N32   
$RUNFDS $INDIR moody_dpdx=-100_N8    
$RUNFDS $INDIR moody_dpdx=-1_N16     
$RUNFDS $INDIR moody_dpdx=-1_N32     
$RUNFDS $INDIR moody_dpdx=-1_N8      
$RUNFDS $INDIR poiseuille_N16_mu025  
$RUNFDS $INDIR poiseuille_N32_mu025  
$RUNFDS $INDIR poiseuille_N64_mu0125 
$RUNFDS $INDIR poiseuille_N64_mu025  
$RUNFDS $INDIR poiseuille_N8_mu025   

$RUNFDS $INDIR z0=p0001_dpdx=-1_N8     
$RUNFDS $INDIR z0=p001_dpdx=-1_N8      
$RUNFDS $INDIR z0=p01_dpdx=-1_N8       
$RUNFDS $INDIR z0=p1_dpdx=-1_N8        
$RUNFDS $INDIR z0=p0001_dpdx=-1_N16    
$RUNFDS $INDIR z0=p001_dpdx=-1_N16     
$RUNFDS $INDIR z0=p01_dpdx=-1_N16      
$RUNFDS $INDIR z0=p0001_dpdx=-100_N8   
$RUNFDS $INDIR z0=p001_dpdx=-100_N8    
$RUNFDS $INDIR z0=p01_dpdx=-100_N8     
$RUNFDS $INDIR z0=p1_dpdx=-100_N8      
$RUNFDS $INDIR z0=p0001_dpdx=-100_N16  
$RUNFDS $INDIR z0=p001_dpdx=-100_N16   
$RUNFDS $INDIR z0=p01_dpdx=-100_N16    
$RUNFDS $INDIR z0=p0001_dpdx=-p0001_N8 
$RUNFDS $INDIR z0=p0001_dpdx=-p01_N8   
$RUNFDS $INDIR z0=p001_dpdx=-p0001_N8  
$RUNFDS $INDIR z0=p001_dpdx=-p01_N8    
$RUNFDS $INDIR z0=p01_dpdx=-p0001_N8   
$RUNFDS $INDIR z0=p01_dpdx=-p01_N8     
$RUNFDS $INDIR z0=p1_dpdx=-p0001_N8    
$RUNFDS $INDIR z0=p1_dpdx=-p01_N8      
$RUNFDS $INDIR z0=p02_dpdx=-1_N16      
$RUNFDS $INDIR z0=p02_dpdx=-100_N16    

echo FDS cases submitted
