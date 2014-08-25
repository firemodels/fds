#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSSC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FSFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_FOSC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SSC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13SFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_13OSC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SSC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24SFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OFW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OFD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OFC.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OSW.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OSD.fds 
$QFDS  $QUEUE -d $INDIR Vettori_Sloped_24OSC.fds 

echo FDS cases submitted
