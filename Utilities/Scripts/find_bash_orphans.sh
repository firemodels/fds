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

cd $SVNROOT/bot
BOTROOT=`pwd`

cd $SVNROOT/smv
SMVROOT=`pwd`

grep $script $FDSROOT/Build/*.html $FDSROOT/Build/*/*.sh $FDSROOT/Manuals/*/*.sh $FDSROOT/Utilities/Scripts/*.sh $FDSROOT/Verification/scripts/*.sh
grep $script $BOTROOT/Bundle/fds/scripts/*.sh $BOTROOT/Bundle/smv/scripts/*.sh $BOTROOT/Firebot/*.sh $BOTROOT/Smokebot/*.sh $BOTROOT/Scripts/*.sh
grep $script $SMVROOT/Build/build*.html $SMVROOT/scripts/*.sh $SMVROOT/Build/*/*.sh $SMVROOT/Utilities/Scripts/*.sh $SMVROOT/Verification/scripts/*.sh
