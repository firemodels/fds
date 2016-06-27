#!/bin/bash
while getopts 'd:' OPTION
do
case $OPTION  in
  d)
   DIR="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

case=$1
cd $DIR
