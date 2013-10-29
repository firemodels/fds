#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`
cd $CURDIR/..

QSMV=/usr/local/bin/qsmokeview.sh
# uncomment following line to stop all cases
#export STOPFDS=1

#rm -f Voltest/frames/mplume8n*.png
#rm -f Voltest/frames/mplumeB8n*.png
rm -f Voltest/frames/voltest2*.png
#rm -f Voltest/frames/voltest3*.png

#$QSMV -d Voltest -p 20 mplume8n
#$QSMV -d Voltest -p 30 -q fire70s  mplumeB8n
$QSMV -d Voltest -p 30 -q vis voltest2
#$QSMV -d Voltest -p 30 voltest3
