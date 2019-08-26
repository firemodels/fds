#!/bin/bash

FED=
MOVIE=
RUNSCRIPT=
ssffile=
dummy=

while getopts 'Ad:fl:mt' OPTION
do
case $OPTION in
  A) # passthrough option
   ;;
  d)
   dir="$OPTARG"
   ;;
  f)
   FED="-fed"
   ;;
  l)
   dummy="$OPTARG"
   ;;
  m)
   MOVIE="y"
   ;;
  t)
   dummy=1
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
echo "--- generating images for: $in.smv, `date`"

scriptfile=$scratchdir/script.$$

notfound=`$SMV -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ];  then
  echo "*** Error: The program $SMV is not available. Run aborted."
  exit
fi
if ! [ -d $fulldir ]; then
  echo "*** Error: The directory $fulldir does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$in.smv ]; then
  echo "*** Error: The smokeview file, $fulldir/$in.smv, does not exist. Run aborted."
  exit
fi
if ! [ -e $fulldir/$ssffile ]; then
  echo "*** Error: The smokeview script file, $fulldir/$ssffile, does not exist. Run aborted."
  exit
fi

#source ~/.bashrc_fds default
cd $fulldir
echo $SMV $FED $SMVBINDIR $RUNSCRIPT $in
$SMV $FED $SMVBINDIR -redirect $RUNSCRIPT $in
