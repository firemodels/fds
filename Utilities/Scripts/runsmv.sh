#!/bin/bash -f

dir=$1
in=$2

fulldir=$BASEDIR/$dir

scriptfile=$scratchdir/script.$$
if ! [ -e $SMV ];  then
  echo "The file $SMV does not exist. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  echo "The smokeview file, $fulldir/$in.smv, does not exist. Run aborted."
  exit
fi

source ~/.bashrc_fds intel64
cd $fulldir
$SMV $SMVBINDIR -runscript $in
