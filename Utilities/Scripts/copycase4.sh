#!/bin/bash

dir=$3
in=$4


cd $dir
if ! [ -d $OUTDIR/$dir ]; then
  mkdir $OUTDIR/$dir
fi

cp $in.fds $OUTDIR/$dir/.
if [ -e $in.ini ]; then
  cp $in.ini $OUTDIR/$dir/.
fi
if [ -e $in.ssf ]; then
  cp $in.ssf $OUTDIR/$dir/.
fi
