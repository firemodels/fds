#!/bin/bash
CURDIR=`pwd`
QUEUE=
INTEL=

function usage {
echo "-h - display this message"
echo "-q queue - specify queue (default: none on any computers "
echo "           without a queuing system)"
echo "-J - use Intel version of fds"
exit
}
notfound=notfound.$$
qstat >& $notfound
status=`cat $notfound | tail -1 | grep "not found" | wc -l`
if [ $status -eq 1 ] ;then
  QUEUE=none
fi

while getopts 'hq:J' OPTION
do
case $OPTION  in
  h)
   usage
   ;;
  J)
   INTEL="-J"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

if [ "$QUEUE" != "" ]; then
  QUEUE="-q $QUEUE"
fi

cd ../../../fds/Verification/scripts
./Run_FDS_Cases.sh $INTEL -j GM_ -g  -W $QUEUE

cd $CURDIR
cd ../../../fds/Verification/scripts
./Make_FDS_Pictures.sh -g
cd $CURDIR

