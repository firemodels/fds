#!/bin/bash

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_100_RI=5.fds
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_150_RI=5.fds
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_200_RI=5.fds
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_300_RI=5.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_400_RI=5.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_500_RI=5.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_50_RI=5.fds  
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_600_RI=5.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_75_RI=5.fds  

$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_100_RI=10.fds
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_150_RI=10.fds
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_200_RI=10.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_300_RI=10.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_400_RI=10.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_500_RI=10.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_50_RI=10.fds  
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_600_RI=10.fds 
$QFDS -p 2 -r $qq -d $INDIR NIST_RSE_1994_75_RI=10.fds  

$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_100_RI=20.fds
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_150_RI=20.fds
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_200_RI=20.fds 
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_300_RI=20.fds 
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_400_RI=20.fds 
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_500_RI=20.fds 
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_50_RI=20.fds  
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_600_RI=20.fds 
$QFDS -p 8 -r $qq -d $INDIR NIST_RSE_1994_75_RI=20.fds
 
echo FDS cases submitted
