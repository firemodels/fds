#!/bin/bash
touch smokebot_running
./smokebot_linux.sh $*
rm smokebot_running
