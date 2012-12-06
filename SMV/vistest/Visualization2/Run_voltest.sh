#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR/..

QFDS=/usr/local/bin/qfds.sh
# uncomment following line to stop all cases
#export STOPFDS=1

$QFDS -r -d Visualization2 plume8c.fds
$QFDS -r -d Visualization2 plume8n.fds
$QFDS -r -d Visualization2 plumeB8c.fds
$QFDS -r -d Visualization2 plumeB8n.fds
$QFDS -r -p 8 -d Visualization2 mplume8c.fds
$QFDS -r -p 8 -d Visualization2 mplume8n.fds
$QFDS -r -p 8 -d Visualization2 mplumeB8c.fds
$QFDS -r -p 8 -d Visualization2 mplumeB8n.fds
