#!/bin/bash -f
dir=$1
in=$2

fulldir=$BASEDIR/$dir

scriptfile=$scratchdir/script.$$
if ! [ -e $SMV ];  then
  echo "The file $SMV does not exit. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  echo "The smokeview input file, $fulldir/$in.smv, does not exit. Run aborted."
  exit
fi

disp=`id -u`
Xvfb :$disp &
rxs=$!

#sleep 8

export DISPLAY=:$disp
source ~/.bashrc_fds intel64
cd $fulldir
$SMV -runscript $in

#sleep 8

kill $rxs
