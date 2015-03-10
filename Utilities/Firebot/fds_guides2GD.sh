#!/bin/bash

CURDIR=`pwd`
PARENT_ID=0B-W-dkXwdHWNeU41UkZiNV93Mmc
FROMDIR=~firebot/firebot/output
BASEDIR=Newest_FDS_Guides
UPDATE_GUIDES=$FROMDIR/$BASEDIR/update_guides

UPLOAD ()
{
  FILE=$1
  FILEnew=${FILE}_new.pdf
  cp $FILE.pdf /tmp/$FILEnew
  drive list  | grep $FILEnew | awk '{ system("drive delete -i " $1)} '
  drive upload -p $PARENT_ID -f /tmp/$FILEnew
}

if [ -e $UPDATE_GUIDES ] ; then
  cd $FROMDIR/$BASEDIR
  rm -f $UPDATE_GUIDES
  UPLOAD FDS_Configuration_Management_Plan
  UPLOAD FDS_Technical_Reference_Guide
  UPLOAD FDS_User_Guide
  UPLOAD FDS_Validation_Guide
  UPLOAD FDS_Verification_Guide
  cd $CURDIR
fi
