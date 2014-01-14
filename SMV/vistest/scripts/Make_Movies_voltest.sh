#!/bin/bash

makemovie=/usr/local/bin/make_movie.sh

CURDIR=`pwd`
export SVNROOT=`pwd`/../../..
cd $SVNROOT
export SVNROOT=`pwd`

OUTDIR=$SVNROOT/Manuals/SMV_Summary/movies

cd $CURDIR/..

#$makemovie -i Voltest/frames -o $OUTDIR mplume8n
#$makemovie -i Voltest/frames -o $OUTDIR mplumeB8n
$makemovie -i Voltest/frames -o $OUTDIR voltest2
