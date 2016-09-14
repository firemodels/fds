#!/bin/bash
script=$1

# This script compiles all LaTeX guides, then checks for "orphaned" files.
# These are files that are not referenced in the LaTeX guides.

CURDIR=`pwd`
SVNROOT=`pwd`/../../..
cd $SVNROOT
SVNROOT=`pwd`

cd $SVNROOT/fds
FDSROOT=`pwd`

cd $SVNROOT/smv
SMVROOT=`pwd`

echo looking for $script in FDS repo
grep $script $FDSROOT/Build/*/*.sh $FDSROOT/Build/Scripts/*.sh $FDSROOT/bot/Firebot/*.sh $FDSROOT/Utilities/Scripts/*.sh $FDSROOT/Verification/scripts/*.sh
echo looking for $script in SMV repo
grep $script $SMVROOT/Build/*/*.sh $SMVROOT/Build/scripts/*.sh $SMVROOT/bot/Smokebot/*.sh $SMVROOT/Utilities/Scripts/*.sh $SMVROOT/Verification/scripts/*.sh
