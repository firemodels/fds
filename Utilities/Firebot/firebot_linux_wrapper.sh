#!/bin/bash

running=~/firebot/firebot_running
if [ -e $running ] ; then
  echo Firebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi
touch $running
svn update
~/firebot/firebot_linux.sh "$@"
rm $running

