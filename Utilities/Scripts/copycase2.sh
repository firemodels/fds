#!/bin/bash -f

dir=$2
in=$3


if ! [ -d $OUTDIR/$dir ]; then
  mkdir $OUTDIR/$dir
fi
cd $dir
cp $in.fds $OUTDIR/$dir/.
if [ -e $in.ini ]; then
  cp $in.ini $OUTDIR/$dir/.
fi
if [ -e $in.ssf ]; then
  cp $in.ssf $OUTDIR/$dir/.
fi
