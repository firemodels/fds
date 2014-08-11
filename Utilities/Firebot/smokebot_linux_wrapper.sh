#!/bin/bash
if [ -e ~/SMOKEBOT/smokebot_running ] ; then
  exit
fi
touch ~/SMOKEBOT/smokebot_running
~/SMOKEBOT/smokebot_linux.sh "$@"
rm ~/SMOKEBOT/smokebot_running
