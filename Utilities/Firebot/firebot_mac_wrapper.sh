#!/bin/bash

DEBUG=
while getopts 'd' OPTION
do
case $OPTION  in
  d)
  DEBUG=-d
  ;;
esac
done
shift $(($OPTIND-1))

running=~/firebot/firebot_running
if [ -e $running ] ; then
  echo A bot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi
svn update
touch $running
~/firebot/firebot_mac_git.sh $DEBUG "$@"
rm $running
