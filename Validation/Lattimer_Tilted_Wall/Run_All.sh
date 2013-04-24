#!/bin/bash -f

export SVNROOT=`pwd`/../..
export QFDS=/usr/local/bin/qfds.sh
export BASEDIR=`pwd`
export INDIR=Current_Results
# qq="-q fire80s"
qq=
source ~/.bashrc_fds intel64

# uncomment following line to stop all cases
# export STOPFDS=1

/bin/sh -c "cp $BASEDIR/FDS_Input_Files/*.fds $BASEDIR/$INDIR"

$QFDS -r -p 27 $qq -d $INDIR Lattimer_20_kW_0_degree_coarse.fds
$QFDS -r -p 27 $qq -d $INDIR Lattimer_20_kW_0_degree.fds
$QFDS -r -p 27 $qq -d $INDIR Lattimer_20_kW_0_degree_fine.fds

$QFDS -r $qq -d $INDIR Lattimer_20_kW_0_degree_ibm.fds
$QFDS -r $qq -d $INDIR Lattimer_20_kW_10_degree_ibm.fds
$QFDS -r $qq -d $INDIR Lattimer_20_kW_20_degree_ibm.fds

echo FDS cases submitted
