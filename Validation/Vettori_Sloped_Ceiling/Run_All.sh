#!/bin/bash -f

export SVNROOT=`pwd`/../..
export FDS=$SVNROOT/FDS_Compilation/intel_linux_64/fds_intel_linux_64
export RUNFDS=$SVNROOT/Utilities/Scripts/runfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$RUNFDS $INDIR Vettori_Sloped_FSSW 
$RUNFDS $INDIR Vettori_Sloped_FSSD 
$RUNFDS $INDIR Vettori_Sloped_FSSC 
$RUNFDS $INDIR Vettori_Sloped_FSFW 
$RUNFDS $INDIR Vettori_Sloped_FSFD 
$RUNFDS $INDIR Vettori_Sloped_FSFC 
$RUNFDS $INDIR Vettori_Sloped_FOFW 
$RUNFDS $INDIR Vettori_Sloped_FOFD 
$RUNFDS $INDIR Vettori_Sloped_FOFC 
$RUNFDS $INDIR Vettori_Sloped_FOSW 
$RUNFDS $INDIR Vettori_Sloped_FOSD 
$RUNFDS $INDIR Vettori_Sloped_FOSC 
$RUNFDS $INDIR Vettori_Sloped_13SSW 
$RUNFDS $INDIR Vettori_Sloped_13SSD 
$RUNFDS $INDIR Vettori_Sloped_13SSC 
$RUNFDS $INDIR Vettori_Sloped_13SFW 
$RUNFDS $INDIR Vettori_Sloped_13SFD 
$RUNFDS $INDIR Vettori_Sloped_13SFC 
$RUNFDS $INDIR Vettori_Sloped_13OFW 
$RUNFDS $INDIR Vettori_Sloped_13OFD 
$RUNFDS $INDIR Vettori_Sloped_13OFC 
$RUNFDS $INDIR Vettori_Sloped_13OSW 
$RUNFDS $INDIR Vettori_Sloped_13OSD 
$RUNFDS $INDIR Vettori_Sloped_13OSC 
$RUNFDS $INDIR Vettori_Sloped_24SSW 
$RUNFDS $INDIR Vettori_Sloped_24SSD 
$RUNFDS $INDIR Vettori_Sloped_24SSC 
$RUNFDS $INDIR Vettori_Sloped_24SFW 
$RUNFDS $INDIR Vettori_Sloped_24SFD 
$RUNFDS $INDIR Vettori_Sloped_24SFC 
$RUNFDS $INDIR Vettori_Sloped_24OFW 
$RUNFDS $INDIR Vettori_Sloped_24OFD 
$RUNFDS $INDIR Vettori_Sloped_24OFC 
$RUNFDS $INDIR Vettori_Sloped_24OSW 
$RUNFDS $INDIR Vettori_Sloped_24OSD 
$RUNFDS $INDIR Vettori_Sloped_24OSC 

echo FDS cases submitted
