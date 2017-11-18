#!/bin/bash

ssffile=
SMOKEVIEW=smokeview
RUNSCRIPT=-runscript
SMOKEVIEWDIR=$(dirname "$0")

while getopts 'd:e:s:' OPTION
do
case $OPTION in
  d)
   dir="$OPTARG"
   ;;
  e)
   SMOKEVIEW="$OPTARG"
   ;;
  s)
   RUNSCRIPT="-scriptfile $OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

in=$1
in=${in%*.*}

ssffile=$in.ssf

notfound=`$SMOKEVIEW -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ];  then
  echo "*** Error (fatal): The program $SMOKEVIEW is not available. Run aborted."
  exit
fi

if ! [ -e $in.smv ]; then
  echo "*** Error (fatal): The smokeview file, $in.smv, does not exist. Run aborted."
  exit
fi
if ! [ -e $ssffile ]; then
  echo "*** Error (fatal): The smokeview script file, $ssffile, does not exist. Run aborted."
  exit
fi


echo $SMOKEVIEW $RUNSCRIPT $in
source $SMOKEVIEWDIR/startXserver.sh
$SMOKEVIEW $RUNSCRIPT $in
source $SMOKEVIEWDIR/stopXserver.sh
