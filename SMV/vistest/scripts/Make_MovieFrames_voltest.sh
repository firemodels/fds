#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR/..

QSMV=/usr/local/bin/qsmokeview.sh
# uncomment following line to stop all cases
#export STOPFDS=1

rm -f Voltest/frames/*.png

#$QSMV -d Voltest -p 30 mplume8n
$QSMV -d Voltest -p 30 mplumeB8n
