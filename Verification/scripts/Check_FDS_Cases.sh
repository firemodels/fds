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
outfile=$DIR/${case%.*}.out
if [  -e $outfile ]; then
  if [[ grep -rI 'completed successfully' $outfile` == "" ]]; then
    # Continue along
  else
    error "*** error: the case $infile did not complete"
  fi
else
  error "*** error: the case $infile did not run ($outfile does not exist)"
fi
