#!/bin/bash
FORCE=$1

gitwebrepo=~/FDS-SMVwebpages
oldpage=~/.firebot/oldpage
newpage=~/.firebot/newpage
olddata=~/.firebot/old_data
newdata=~/.firebot/fds_times.csv
running=~/.fdssmvgit/bot_running
curdir=`pwd`
EXIT="yes"

# don't update status page if firebot is running
if [ -e $running ] ; then
  if [ "$FORCE" == "" ]; then
    exit
  fi
fi

# check if status web page has changed

./make_pubpage.sh -b > $newpage
if [ ! -e $oldpage ]; then
  cp $newpage $oldpage
fi
ndiff=`diff $oldpage $newpage|wc -l`
if [ ! "$ndiff" == "0" ] ; then
  cp $newpage $oldpage
  EXIT="no"
fi

# check if FDS benchmark times have changed

./make_timelist.sh > $newdata
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

./make_pubpage.sh > $newpage

cd $gitwebrepo
git remote update
git merge origin/nist-pages

cp $newpage firebot_status.html
git add firebot_status.html

git commit -m "firebot: update firebot status page `date`"
git push
