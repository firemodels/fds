#!/bin/bash
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
  exit
fi

# check if status web page has changed

./smvstatus2html.sh  > $newpage
if [ ! -e $oldpage ]; then
  cp $newpage $oldpage
fi
ndiff=`diff $oldpage $newpage|wc -l`
if [ ! "$ndiff" == "0" ] ; then
  cp $newpage $oldpage
  EXIT="no"
fi

# check if FDS benchmark times have changed

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
   exit
fi

./smvstatus_pubtop.sh > $newpage
./smvstatus2html.sh >> $newpage
./status_pubbot.sh >> $newpage

cd $gitwebrepo
git remote update
git merge origin/gh-pages

cp $newpage smokebot_status.html
git add smokebot_status.html

cp $newdata smv_times.csv
git add smv_times.csv

git commit -m "smokebot: update smokebot status page `date`"
git push
