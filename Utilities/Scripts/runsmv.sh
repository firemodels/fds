#!/bin/bash -f

function smv_linux {
disp=`id -u`
Xvfb :$disp -fp /usr/share/X11/fonts/misc -screen 0 1280x1024x24 &
rxs=$!
echo pausing for 8 seconds
sleep 8

export DISPLAY=:$disp
source ~/.bashrc_fds intel64
cd $fulldir
$SMV $SMVBINDIR -runscript $in
echo pausing for 8 seconds
sleep 8

kill $rxs
}

function smv_osx {
source ~/.bashrc_fds intel64
cd $fulldir
$SMV $SMVBINDIR -runscript $in
}

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

OS=`uname`
if [ "$OS" == "Darwin" ]; then
  OSSYSTEM=osx
else
  OSSYSTEM=linux
fi

if [ "$OSSYSTEM" == "linux" ]; then
  smv_linux;
else
  smv_osx;
fi
