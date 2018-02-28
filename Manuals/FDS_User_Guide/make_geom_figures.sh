#!/bin/bash
CURDIR=`pwd`
QUEUE=

function usage {
echo "-h - display this message"
echo "-q queue - specify queue"
exit
}

while getopts 'hq:' OPTION
do
case $OPTION  in
  h)
   usage
   ;;
  q)
   QUEUE="-q $OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

cd ../../../fds/Verification/scripts
#./Run_FDS_Cases.sh -J -j GM_ -g  -W $QUEUE

cd $CURDIR
cd ../../../fds/Verification/scripts
./Make_FDS_Pictures.sh -g
cd $CURDIR

