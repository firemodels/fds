#!/bin/bash
SVNROOT=~/FDS-SMV
mailTo="gforney@gmail.com koverholt@gmail.com"

SVNSOURCE=$SVNROOT/SMV/source
SVNFILE=$SVNROOT/svnversion
SMOKEBOTDIR=~/SMOKEBOT/
SMOKEBOTEXE=./run_smokebot.sh

cd $SVNSOURCE
svn update
THISSVN=`svn info | tail -3 | head -1 | awk '{print $4}'`
LASTSVN=`cat $SVNFILE`
if [ $THISSVN != $LASTSVN ] ; then
  echo $THISSVN>$SVNFILE
  echo -e "smokeview source has changed\n (Old revision: $LASTSVN, New revision: $THISSVN)\n Smokebot run initiated.\n" | mail -s "smokebot run initiated" $mailTo > /dev/null
  cd $SMOKEBOTDIR
  $SMOKEBOTEXE
fi
