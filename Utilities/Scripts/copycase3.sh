#!/bin/bash -f

dir=$1
in=$2


cd $dir
if ! [ -d $OUTDIR/$dir ]; then
  mkdir $OUTDIR/$dir
fi

cp $in.in $OUTDIR/$dir/.
if [ -e $in.ini ]; then
  cp $in.ini $OUTDIR/$dir/.
fi
if [ -e $in.ssf ]; then
  cp $in.ssf $OUTDIR/$dir/.
fi
