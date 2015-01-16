#!/bin/bash

DEBUG=
REVERT=
while getopts 'dn' OPTION
do
case $OPTION  in
  d)
  $DEBUG=-d
  ;;
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
~/firebot/firebot_mac.sh $DEBUG $REVERT "$@"
rm $running
