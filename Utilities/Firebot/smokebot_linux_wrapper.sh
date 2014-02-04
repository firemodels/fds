#!/bin/bash
touch ~/SMOKEBOT/smokebot_running
~/SMOKEBOT/smokebot_linux.sh $*
rm ~/SMOKEBOT/smokebot_running
