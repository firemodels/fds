#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=

# uncomment following line to stop all cases
#export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r $qq -d $INDIR NIST_He_14400_LC_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_LC_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_LC_ULV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_LR_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_LR_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_LR_ULV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_UC_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_UC_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_14400_UC_ULV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LC_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LC_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LC_ULV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LR_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LR_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_LR_ULV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_UC_SLV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_UC_SSV.fds
$QFDS -r $qq -d $INDIR NIST_He_3600_UC_ULV.fds

echo FDS cases submitted
