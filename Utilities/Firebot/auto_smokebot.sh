#!/bin/bash
SVNROOT=~/FDS-SMV
mailTo="gforney@gmail.com koverholt@gmail.com"

SMV_SOURCE=$SVNROOT/SMV/source
SVN_SMVFILE=$SVNROOT/smv_revision

FDS_SOURCE=$SVNROOT/FDS_Source
SVN_FDSFILE=$SVNROOT/fds_revision

SMOKEBOTDIR=~/SMOKEBOT/
SMOKEBOTEXE=./run_smokebot.sh

MESSAGE_FILE=$SVNROOT/message

cd $SMV_SOURCE
svn update > /dev/null
THIS_SMVSVN=`svn info | tail -3 | head -1 | awk '{print $4}'`
THIS_SMVAUTHOR=`svn info | tail -4 | head -1 | awk '{print $4}'`
LAST_SMVSVN=`cat $SVN_SMVFILE`

cd $FDS_SOURCE
svn update > /dev/null
THIS_FDSSVN=`svn info | tail -3 | head -1 | awk '{print $4}'`
THIS_FDSAUTHOR=`svn info | tail -4 | head -1 | awk '{print $4}'`
LAST_FDSSVN=`cat $SVN_FDSFILE`

rm -f $MESSAGE_FILE
if [[ $THIS_SMVSVN != $LAST_SMVSVN || $THIS_FDSSVN != $LAST_FDSSVN ]] ; then
  if [[ $THIS_SMVSVN != $LAST_SMVSVN ]] ; then
    echo -e "smokeview source has changed. $LAST_SMVSVN->$THIS_SMVSVN($THIS_SMVAUTHOR)" >> $MESSAGE_FILE
  fi
  if [[ $THIS_FDSSVN != $LAST_FDSSVN ]] ; then
    echo -e "FDS source has changed. $LAST_FDSSVN->$THIS_FDSSVN($THIS_FDSAUTHOR)" >> $MESSAGE_FILE
  fi
  echo -e "Smokebot run initiated." >> $MESSAGE_FILE
  echo $THIS_SMVSVN>$SVN_SMVFILE
  echo $THIS_FDSSVN>$SVN_FDSFILE
  cat $MESSAGE_FILE | mail -s "smokebot run initiated" $mailTo > /dev/null
  cd $SMOKEBOTDIR
  $SMOKEBOTEXE
fi
