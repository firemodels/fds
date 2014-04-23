#!/bin/bash
# $Date: 2012-11-28 11:04:24 -0500 (Wed, 28 Nov 2012) $ 
# $Revision: 13952 $
# $Author: gforney $

PROG=$0

# setup default queue name

progname=qsmokeview.sh

nprocesses=1
dir=.
queue=fire70s

if [ $# -lt 1 ]
then
  echo "Usage: $progname [-p nprocesses] casename "
  echo ""
  echo " -d directory - directory containing case"
  echo " -q queue - queue used to run cases (only fire70s or vis)"
  echo " -p nprocesses - number of processes used to run a case [default: 1] "
  echo "input_file - input file"
  echo ""
  exit
fi

VOLRENDER="-x "
SCRIPTFILE=

# default parameter settings

while getopts 'd:m:p:q:w' OPTION
do
case $OPTION  in
  d)
  dir="$OPTARG"
  ;;
  m)
   SCRIPTFILE="$OPTARG"
  ;;
  p)
   nprocesses="$OPTARG"
  ;;
  q)
  queue="$OPTARG"
  ;;
  w)
  VOLRENDER=
  ;;
esac
done
shift $(($OPTIND-1))
input=$1

if [ "$SCRIPTFILE" != "" ] ; then
  SCRIPTFILE="-m $SCRIPTFILE"
fi
echo SCRIPTFILE=$SCRIPTFILE

stopfile=$dir/$input.stop

if [ $STOPFDS ]; then
 echo "stopping case: $input"
 touch $stopfile
 exit
fi

if [ "$queue" == "batch" ] ; then
  queue=fire70s
fi

i=0
while [ $i -lt $nprocesses ]; do
  qfds.sh -d $dir -r -s $VOLRENDER $SCRIPTFILE -y $i -z $nprocesses -q $queue $input
  sleep 5
  let i=i+1
done
