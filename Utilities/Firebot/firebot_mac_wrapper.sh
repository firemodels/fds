#!/bin/bash

DEBUG=
while getopts 'dn' OPTION
do
case $OPTION  in
  d)
  $DEBUG=-d
  ;;
esac
done
shift $(($OPTIND-1))

running=~/firebot/firebot_running
if [ -e $running ] ; then
  exit
fi
touch $running
~/firebot/firebot_mac.sh $DEBUG "$@"
rm $running
