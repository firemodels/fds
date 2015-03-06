#!/bin/bash

RUNSCRIPT=
ssffile=
SMV=smokeview
SMVDIR=$(dirname "$0")

while getopts 'd:' OPTION
do
case $OPTION in
  d)
   dir="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

in=$1
in=${in%*.*}

RUNSCRIPT=-runscript
ssffile=$in.ssf

notfound=`$SMV -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ];  then
  echo "*** Error (fatal): The program $SMV is not available. Run aborted."
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


echo $SMV $RUNSCRIPT $in
source $SMVDIR/startXserver.sh
$SMV $RUNSCRIPT $in
source $SMVDIR/stopXserver.sh
