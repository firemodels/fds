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

$QFDS -r -p 16 $qq -d $INDIR Sandia_He_1m_dx6cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_He_1m_dx3cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_He_1m_dx1p5cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test14_dx6cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test14_dx3cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test14_dx1p5cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test17_dx6cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test17_dx3cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test17_dx1p5cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test24_dx6cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test24_dx3cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_CH4_1m_Test24_dx1p5cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_H2_1m_Test35_dx6cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_H2_1m_Test35_dx3cm.fds
$QFDS -r -p 16 $qq -d $INDIR Sandia_H2_1m_Test35_dx1p5cm.fds

echo FDS cases submitted
