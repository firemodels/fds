#!/bin/bash
CURDIR=`pwd`
FROMDIR=~smokebot/smokebot/output
BASEDIR=Newest_Smokeview_Guides
UPDATE_GUIDES=$FROMDIR/$BASEDIR/update_guides
if [ -e $UPDATE_GUIDES ] ; then
  cd $FROMDIR
  gdrive upload -i $BASEDIR
  rm $UPDATE_GUIDES
  cd $CURDIR
fi
