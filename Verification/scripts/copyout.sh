#!/bin/bash

outdir=$FDS_SVNROOT/Verification/scripts/Outfiles
dir=$1
infile=$2

fulldir=$FDS_SVNROOT/Verification/$dir
outfile=$infile.out

# ensure that various files and directories exist

if ! [ -d $fulldir ]; then
  exit
fi
if ! [ -e $fulldir/$outfile ]; then
  exit
fi
ln -f -s $fulldir/$outfile $outdir/$outfile
