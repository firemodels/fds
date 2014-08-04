#!/bin/bash

running=~/firebot/firebot_running
if [ -e $running ] ; then
  exit
fi
touch $running
~/firebot/firebot_mac.sh $*
rm $running
