#!/bin/bash -f

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/..
export SMV=~/bin/smv5_linux_64
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export BASEDIR=`pwd`

./FDS_Picture_Cases.sh

echo FDS cases submitted

