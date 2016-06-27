#!/bin/bash
DIR=
while getopts 'd:tp:o:' OPTION
do
case $OPTION in
  d)
  DIR="$OPTARG"
  ;;
  o)
  dummy="$OPTARG"
  ;;
  p)
  dummy2="$OPTARG"
  ;;
  t)
  ;;
esac
done
shift $(($OPTIND-1))

CASE=${1%.*}.out
CPUTIME=`grep Total ../$DIR/$CASE 2>/dev/null| grep Elapsed | grep Wall | tail -1 | awk '{print $7}'`
echo $CPUTIME $DIR/$CASE 

