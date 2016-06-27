#!/bin/bash

# parse options
while getopts 'Abd:e:im:n:o:p:q:st' OPTION
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
rm -f $case*.s3d
rm -f $case*.sf
rm -f $case*.bf
rm -f $case*.iso
rm -f $case*.sz
rm -f $case*.csv
rm -f $case*.q
rm -f $case*.restart
rm -f $case*.out
rm -f $case*.err
