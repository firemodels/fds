#!/bin/bash
script=$1

CURDIR=`pwd`
SVNROOT=`pwd`/../../..
cd $SVNROOT
SVNROOT=`pwd`

cd $SVNROOT/fds
FDSROOT=`pwd`

cd $SVNROOT/smv
SMVROOT=`pwd`

echo looking for $script in FDS repo
grep $script $FDSROOT/Build/*.html $FDSROOT/Utilities/*.html $FDSROOT/Build/*/*.bat $FDSROOT/Build/Scripts/*.bat $FDSROOT/bot/Firebot/*.bat $FDSROOT/Utilities/Scripts/*.bat $FDSROOT/Verification/scripts/*.bat

echo looking for $script in SMV repo
grep $script $SMVROOT/Build/*.html $SMVROOT/Build/*/*.bat $SMVROOT/Build/Scripts/*.bat $SMVROOT/bot/Firebot/*.bat $SMVROOT/Utilities/Scripts/*.bat $SMVROOT/Verification/scripts/*.bat
