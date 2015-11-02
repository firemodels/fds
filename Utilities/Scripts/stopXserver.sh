#!/bin/bash
SSH=
while getopts 'S:' OPTION
do
case $OPTION  in
  S)
   SSH="ssh $OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

echo shutting down graphics environment
sleep 8
if [ "`uname`" == "Darwin" ]; then
  PIDS=`$SSH ps -u $USER | grep Xvfb | grep -v grep |  awk '{print $2}'`
  for f in $PIDS
  do
    $SSH kill -9 $f
  done
else
  $SSH kill $SMV_ID
fi

