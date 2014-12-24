#!/bin/bash

# parse options
DIR=.
while getopts 'bd:e:im:n:o:p:q:st' OPTION
do
case $OPTION  in
  b)
   dummy=1
   ;;
  d)
   DIR="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))


outdir=$FDS_SVNROOT/Verification/scripts/Outfiles
in=$1
infile=${in%.*}

fulldir=$FDS_SVNROOT/Verification/$DIR
outfile=$infile.out

# ensure that various files and directories exist

if ! [ -d $fulldir ]; then
  exit
fi
if ! [ -e $fulldir/$outfile ]; then
  exit
fi
ln -f -s $fulldir/$outfile $outdir/$outfile
