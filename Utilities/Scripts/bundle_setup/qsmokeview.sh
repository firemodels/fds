#!/bin/bash
# $Date: 2012-11-28 11:04:24 -0500 (Wed, 28 Nov 2012) $ 
# $Revision: 13952 $
# $Author: gforney $

PROG=$0

# setup default queue name

progname=qsmokeview.sh

nprocesses=1

if [ $# -lt 1 ]
then
  echo "Usage: $progname [-p nprocesses] casename "
  echo ""
  echo " -p nprocesses - number of processes used to run a case [default: 1] "
  echo "input_file - input file"
  echo ""
  exit
fi

# default parameter settings

while getopts 'p:' OPTION
do
case $OPTION  in
  p)
   nprocesses="$OPTARG"
  ;;
esac
done
shift $(($OPTIND-1))
input=$1

i=0
while [ $i -lt $nprocesses ]; do
  qfds.sh -r -s -x -y $i -z $nprocesses -q fire70s $input
  sleep 3
  let i=i+1
done
