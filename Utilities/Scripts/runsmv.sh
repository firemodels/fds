#!/bin/bash

FED=
MOVIE=
RUNSCRIPT=
ssffile=
WFDSCASE="no"
TIMEOPTION=

while getopts 'd:fmtw' OPTION
do
case $OPTION in
  d)
   dir="$OPTARG"
   ;;
  f)
   FED="-fed"
   ;;
  m)
   MOVIE="y"
   ;;
  t)
   TIMEOPTION=-time
   ;;
  w)
   WFDSCASE="y"
   ;;
esac
done
shift $(($OPTIND-1))

in=$1
in=${in%*.*}

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

notfound=`$SMV -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ];  then
  echo "*** Error (fatal): The program $SMV is not available. Run aborted."
  exit
fi

if ! [ -d $fulldir ]; then
  echo "*** Error (fatal): The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  if [ "$WFDSCASE" == "n" ]; then
    echo "*** Error (fatal): The smokeview file, $fulldir/$in.smv, does not exist. Run aborted."
  else
    echo "Warning: The smokeview file, $fulldir/$in.smv, does not exist."
  fi
  exit
fi
if ! [ -e $fulldir/$ssffile ]; then
  if [ "$WFDSCASE" == "n" ]; then
    echo "*** Error (fatal): The smokeview script file, $fulldir/$ssffile, does not exist. Run aborted."
  else
    echo "Warning: The smokeview script file, $fulldir/$ssffile, does not exist."
  fi
  exit
fi

source ~/.bashrc_fds
cd $fulldir
echo $SMV $FED $SMVBINDIR $RUNSCRIPT $in
$SMV $TIMEOPTION $FED $SMVBINDIR -redirect $RUNSCRIPT $in
