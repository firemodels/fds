#!/bin/csh -f
echo updating repository
cd ~/FDS-SMV
svn update
echo starting smokebot
cd ~/SMOKEBOT
./run_smokebot.sh -m
