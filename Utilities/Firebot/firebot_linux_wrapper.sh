#!/bin/bash

running=~/firebot/firebot_running
if [ -e $running ] ; then
  exit
fi
touch $running
svn update
~/firebot/firebot_linux.sh "$@"
rm $running

