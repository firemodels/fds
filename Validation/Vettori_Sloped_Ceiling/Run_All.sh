#!/bin/bash

# This script runs a set of Validation Cases on a Linux machine with a batch queuing system.
# See the file Validation/Common_Run_All.sh for more information.
export SVNROOT=`pwd`/../..
source $SVNROOT/Validation/Common_Run_All.sh

$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFC.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSW.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSD.fds 
$QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSC.fds 

# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSSC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FSFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_FOSC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SSC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13SFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_13OSC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SSC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24SFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OFC_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSW_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSD_geom.fds 
# $QFDS $DEBUG $QUEUE -d $INDIR Vettori_Sloped_24OSC_geom.fds 

echo FDS cases submitted
