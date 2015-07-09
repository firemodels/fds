#!/bin/bash
if [ -e smokebot_running ] ; then
  exit
fi
touch smokebot_running
./smokebot_linux_git.sh "$@"
rm ./smokebot_running
