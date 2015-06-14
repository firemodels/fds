#!/bin/bash

firebotdir=~/firebotgit
running=$firebotdir/firebot_running
if [ -e $running ] ; then
  echo Firebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi
touch $running
$firebotdir/firebot_linux_git.sh "$@"
rm $running

