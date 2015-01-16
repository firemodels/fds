#!/bin/bash

REVERT=
while getopts 'n' OPTION
do
case $OPTION  in
  n)
  REVERT=-n
  ;;
esac
done
shift $(($OPTIND-1))

running=~/firebot/firebot_running
if [ -e $running ] ; then
  exit
fi
touch $running
~/firebot/firebot_mac.sh $REVERT "$@"
rm $running
