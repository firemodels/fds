#!/bin/bash
CURDIR=`pwd`
FROMDIR=~firebot/firebot/output
BASEDIR=Newest_FDS_Guides
UPDATE_GUIDES=$FROMDIR/$BASEDIR/update_guides
if [ -e $UPDATE_GUIDES ] ; then
  cd $FROMDIR
  gdrive upload -i $BASEDIR
  rm $UPDATE_GUIDES
  cd $CURDIR
fi
