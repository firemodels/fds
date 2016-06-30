#!/bin/bash
FORCE=$1

gitwebrepo=~/FDS-SMVgitweb
firebotdir=~/FDS-SMVgitclean/Utilities/Firebot
oldpage=~/.smokebot/oldpage
newpage=~/.smokebot/newpage
olddata=~/.smokebot/old_data
newdata=~/.smokebot/smv_times.csv
running=~/.fdssmvgit/bot_running
curdir=`pwd`
EXIT="yes"

# don't update status page if firebot is running
cd $firebotdir
if [ -e $running ] ; then
  if [ "$FORCE" == "" ]; then
    exit
  fi
fi

# check if status web page has changed

./make_pubpage.sh -s -b  > $newpage
if [ ! -e $oldpage ]; then
  cp $newpage $oldpage
fi
ndiff=`diff $oldpage $newpage|wc -l`
if [ ! "$ndiff" == "0" ] ; then
  cp $newpage $oldpage
  EXIT="no"
fi

# check if FDS benchmark times have changed

./make_timelist.sh -s > $newdata
if [ ! -e $olddata ]; then
  cp $newdata $olddata
fi
ndiff=`diff $olddata $newdata|wc -l`
if [ ! "$ndiff" == "0" ] ; then
   cp $newdata $olddata
   EXIT="no"
fi

# if nothing has changed then exit without committing any files
if [ "$EXIT" == "yes" ]; then
  if [ "$FORCE" == "" ]; then
    exit
  fi
fi

./make_pubpage.sh -s > $newpage

cd $gitwebrepo
git remote update
git merge origin/nist-pages

cp $newpage smokebot_status.html
git add smokebot_status.html

git commit -m "smokebot: update smokebot status page `date`"
git push
