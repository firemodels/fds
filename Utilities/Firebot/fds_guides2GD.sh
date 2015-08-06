#!/bin/bash

GDRIVE=~/bin/gdrive
CURDIR=`pwd`
# directory containing guides on google drive : FDS-SMV Newest Manuals
PARENT_ID=0B_wB1pJL2bFQUlJwMmNfaHlqME0

# directory containing guides within firebot account
FROMDIR=~firebot/firebotgit/output
BASEDIR=Newest_Guides

UPLOAD ()
{
  FILE=$1
  FILEnew=${FILE}_new.pdf
  cp $FILE.pdf /tmp/$FILEnew
  $GDRIVE list  | grep $FILEnew | awk '{ system("~/bin/gdrive delete -i " $1)} '
  $GDRIVE upload -p $PARENT_ID -f /tmp/$FILEnew
}

if [ -e $GDRIVE ] ; then
  cd $FROMDIR/$BASEDIR
  UPLOAD FDS_Configuration_Management_Plan
  UPLOAD FDS_Technical_Reference_Guide
  UPLOAD FDS_User_Guide
  UPLOAD FDS_Validation_Guide
  UPLOAD FDS_Verification_Guide
  cd $CURDIR
fi
