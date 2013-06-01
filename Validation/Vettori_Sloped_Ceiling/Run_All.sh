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

$QFDS -r $qq -d $INDIR Vettori_Sloped_FSSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FSSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FSSC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FSFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FSFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FSFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_FOSC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SSC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13SFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_13OSC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SSC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24SFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OFW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OFD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OFC.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OSW.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OSD.fds 
$QFDS -r $qq -d $INDIR Vettori_Sloped_24OSC.fds 

echo FDS cases submitted
