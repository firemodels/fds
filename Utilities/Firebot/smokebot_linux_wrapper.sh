#!/bin/bash
if [ -e smokebot_running ] ; then
  exit
fi
touch smokebot_running
svn update
./smokebot_linux.sh "$@"
rm ./smokebot_running
