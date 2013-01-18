#!/bin/bash
CURDIR=`pwd`
export SVNROOT=`pwd`/../../..

cd $SVNROOT
export SVNROOT=`pwd`

makemovie=/usr/local/bin/make_movie.sh
# uncomment following line to stop all cases
#export STOPFDS=1

cd $CURDIR/..

#$makemovie -i Voltest/frames -o movies mplume8n
#$makemovie -i Voltest/frames -o movies mplumeB8n
$makemovie -i Voltest/frames -o movies voltest2
