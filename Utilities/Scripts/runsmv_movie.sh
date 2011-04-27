#!/bin/bash -f
dir=$1
in=$2
movies=_movies
smvscript=$in$movies.ssf

fulldir=$BASEDIR/$dir

if ! [ -e $SMV ];  then
  echo "The file $SMV does not exit. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "The directory $fulldir does not exit. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  echo "The smokeview file, $fulldir/$in.smv, does not exit. Run aborted."
  exit
fi

disp=`id -u`
Xvfb :$disp -screen 0 1280x1024x24 &
rxs=$!

sleep 8

export DISPLAY=:$disp
source ~/.bashrc_fds intel64
cd $fulldir
echo fulldir=$fulldir
echo 
$SMV -script $smvscript $in

sleep 8

kill $rxs
