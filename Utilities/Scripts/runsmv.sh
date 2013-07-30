#!/bin/bash

FED=
MOVIE=
RUNSCRIPT=
ssffile=

while getopts 'fm' OPTION
do
case $OPTION in
  f)
   FED="-fed"
   ;;
  m)
   MOVIE="y"
   ;;
esac
done
shift $(($OPTIND-1))

dir=$1
in=$2

if [ "$FED" == "" ]; then
if [ "$MOVIE" == "" ]; then
  RUNSCRIPT=-runscript
  ssffile=$in.ssf
else
  MOVIE=_movies
  RUNSCRIPT="-script $in$MOVIE.ssf"
  ssffile=$in$MOVIE.ssf
fi
fi

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
if ! [ -e $fulldir/$ssffile ]; then
  echo "*** Error (fatal): The smokeview script file, $fulldir/$ssffile, does not exist. Run aborted."
  exit
fi

source ~/.bashrc_fds intel64
cd $fulldir
echo $SMV $FED $SMVBINDIR $RUNSCRIPT $in
$SMV $FED $SMVBINDIR -redirect $RUNSCRIPT $in
