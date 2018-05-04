#!/bin/bash

DIR=.

# parse options
while getopts 'ACd:e:Ef:hHiILm:MNn:o:O:p:Pq:rsStT:vVw:' OPTION
do
case $OPTION  in
  d)
   DIR="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

case=$1
infile=$DIR/${case%.*}.fds
errfile=$DIR/${case%.*}.err
if [  -e $errfile ]; then
  if [[ `grep -iE 'completed successfully|Set-up only|stopped by KILL control|TGA analysis only' $errfile` != "" ]]; then
    # continue along
    :
  else
    echo "ERROR: $infile started but did not complete"
  fi
else
  echo "ERROR: $infile did not run"
fi
