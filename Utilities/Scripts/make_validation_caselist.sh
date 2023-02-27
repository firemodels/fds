#!/bin/bash

# This script generates a list of validation test cases to be run
MAKEGITENTRY(){
DIR=$1
gitrevisions=/tmp/gitrevisions.$$

cat $GIT_ROOT/out/$DIR/*git.txt 2> /dev/null | sort -u > $gitrevisions
gitrev=`head -1 $gitrevisions`
if [ "$gitrev" != "" ] ; then
  gitrevshort=`echo $gitrev | awk -F - '{print $3}' | sed 's/^.\{1\}//'`
  gitdate=`git show -s --format=%aD $gitrevshort 2> /dev/null | head -1 | awk '{print $3,$2",",$4}'`
  if [ "$gitdate" == "" ]; then
    gitdate="undefined"
    gitdate2=2000000000
  else
    gitdate2=`git show -s --format=%at $gitrevshort | head -1 | awk '{print $1}'`
  fi
  echo "$DIR ! $gitdate2 ! $gitdate  ! $gitrev  "
fi
rm $gitrevisions
}


CURRENT_DIR=`pwd`
SCRIPTDIR=`dirname "$(readlink -f "$0")"`
cd $SCRIPTDIR/../../..
FIREMODELS_ROOT=`pwd`
export FIREMODELS_ROOT
cd $CURRENT_DIR

while getopts 'r:' OPTION
do
case $OPTION  in
  r)
   FIREMODELS_ROOT="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

cd $FIREMODELS_ROOT/fds/Utilities/Scripts

makelist=$FIREMODELS_ROOT/fds/Validation/Process_All_Output.sh
CASELIST=/tmp/temp.out.$$
TABLE_ENTRIES=/tmp/temp2.out.$$
grep PROCESS $makelist | awk 'BEGIN { FS = " " } ; { print $2 }' | awk '{if(NR>1)print}'> $CASELIST
while read p; do
  MAKEGITENTRY   $p  >> $TABLE_ENTRIES
done <$CASELIST

cat $TABLE_ENTRIES | sort -t"!" -n -k2,2 | awk -F "!" '{ print $1 "!" $2 "!" $3}'

rm $CASELIST $TABLE_ENTRIES


cd $CURRENT_DIR

