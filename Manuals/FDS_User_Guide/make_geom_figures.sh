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

cd ../../Verification/scripts
./Run_SMV_Cases.sh -g  -w -j GM_ $QUEUE

cd $CURDIR
cd ../../Verification/scripts
./Make_SMV_Pictures.sh -g
cd $CURDIR

