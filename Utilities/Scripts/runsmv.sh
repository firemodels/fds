#!/bin/bash -f

dir=$1
in=$2

fulldir=$BASEDIR/$dir
echo ""
echo "--- generating images for: $in.smv"

scriptfile=$scratchdir/script.$$
if ! [ -e $SMV ];  then
  echo "*** Error (fatal): The file $SMV does not exist. Run aborted"
  exit
fi
if ! [ -d $fulldir ]; then
  echo "*** Error (fatal): The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  echo "*** Error (fatal): The smokeview file, $fulldir/$in.smv, does not exist. Run aborted."
  exit
fi

source ~/.bashrc_fds intel64
cd $fulldir
$SMV $SMVBINDIR -redirect -runscript $in
