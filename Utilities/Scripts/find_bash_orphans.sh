#!/bin/bash
script=$1

# This script compiles all LaTeX guides, then checks for "orphaned" files.
# These are files that are not referenced in the LaTeX guides.

CURDIR=`pwd`
SVNROOT=`pwd`/../../..
cd $SVNROOT
SVNROOT=`pwd`

cd $SVNROOT/FDS
FDSROOT=`pwd`

cd $SVNROOT/SMV
SMVROOT=`pwd`

echo looking for $script in FDS repo
grep $script $FDSROOT/Build/*/*.sh $FDSROOT/Build/Scripts/*.sh $FDSROOT/Utilities/Firebot/*.sh $FDSROOT/Utilities/Scripts/*.sh $FDSROOT/Verification/scripts/*.sh
echo looking for $script in SMV repo
grep $script $SMVROOT/Build/*/*.sh $SMVROOT/Build/scripts/*.sh $SMVROOT/Utilities/Smokebot/*.sh $SMVROOT/Utilities/Scripts/*.sh $SMVROOT/Verification/scripts/*.sh
