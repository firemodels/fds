#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR/..

QFDS=/usr/local/bin/qfds.sh
# uncomment following line to stop all cases
export STOPFDS=1

$QFDS -r -d Voltest plume8c.fds
$QFDS -r -d Voltest plume8n.fds
$QFDS -r -d Voltest plumeB8c.fds
$QFDS -r -d Voltest plumeB8n.fds
$QFDS -r -d Voltest voltest1.fds
$QFDS -r -p 8 -d Voltest mplume8c.fds
$QFDS -r -p 8 -d Voltest mplume8n.fds
$QFDS -r -p 8 -d Voltest mplumeB8c.fds
$QFDS -r -p 8 -d Voltest mplumeB8n.fds
$QFDS -r -p 8 -d Voltest mvoltest1.fds
