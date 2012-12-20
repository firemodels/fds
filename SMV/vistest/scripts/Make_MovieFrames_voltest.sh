#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR/..

QSMV=/usr/local/bin/qsmokeview.sh
# uncomment following line to stop all cases
#export STOPFDS=1

rm -f Voltest/frames/mplumeB8n*.png

#$QSMV -d Voltest -p 20 mplume8n
$QSMV -d Voltest -p 50 -q fire7080s  mplumeB8n
#$QSMV -d Voltest -p 30 voltest2
