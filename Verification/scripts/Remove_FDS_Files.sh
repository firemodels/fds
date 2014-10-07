#!/bin/bash

# parse options
while getopts 'bd:e:im:n:o:p:q:st' OPTION
do
case $OPTION  in
  b)
   dummy=1
   ;;
  d)
   DIR="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

case=$1
cd $DIR
rm -f $case*.s3d >& /dev/null
rm -f $case*.sf >& /dev/null
rm -f $case*.bf >& /dev/null
rm -f $case*.iso >& /dev/null
rm -f $case*.sz >& /dev/null
rm -f $case*.csv >& /dev/null
rm -f $case*.q >& /dev/null
rm -f $case*.restart /dev/null
rm -f $case*.out /dev/null
rm -f $case*.err /dev/null
