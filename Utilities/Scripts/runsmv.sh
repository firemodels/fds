#!/bin/bash -f

function smv_linux {
source ~/.bashrc_fds intel64
cd $fulldir
$SMV $SMVBINDIR -runscript $in
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
