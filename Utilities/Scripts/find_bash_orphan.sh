#!/bin/bash
file=$1
type=$2

CURDIR=`pwd`
cd ../../..
GITROOT=`pwd`
BOTREPO=$GITROOT/bot
FDSREPO=$GITROOT/fds
SMVREPO=$GITROOT/smv

echo "************************************************"
echo "*** looking for $file in $BOTREPO"
echo "************************************************"
cd $BOTREPO
grep -r $file --include=*sh

cd $FDSREPO
echo "************************************************"
echo "*** looking for $file in $FDSREPO"
echo "************************************************"
grep -r $file --include=*sh

cd $SMVREPO
echo "************************************************"
echo "*** looking for $file in $SMVREPO"
echo "************************************************"
grep -r $file --include=*sh
