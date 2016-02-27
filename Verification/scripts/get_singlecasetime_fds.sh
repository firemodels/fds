#!/bin/bash
DIR=
while getopts 'Ad:tp:o:' OPTION
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
CPUTIME=`grep MAIN $DIR/$CASE 2>/dev/null| head -1 | awk '{print $2}'`
echo $CPUTIME $DIR/$CASE 

