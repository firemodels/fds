#!/bin/bash -f

# This script runs the FDS Verification Cases on a linux machine with
# a batch queuing system

export SVNROOT=`pwd`/..
export SMV=$SVNROOT/SMV/Build/intel_linux_test_64/smokeview_linux_test_64
export RUNSMV=$SVNROOT/Utilities/Scripts/runsmv.sh
export SMVBINDIR="-bindir $SVNROOT/SMV/for_bundle/"
export BASEDIR=`pwd`

./FDS_Pictures.sh

echo FDS cases submitted

