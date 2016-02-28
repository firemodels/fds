#!/bin/bash
gitwebrepo=~/FDS-SMVgitweb
firebotdir=~/FDS-SMVgitclean/Utilities/Firebot
old=~/.firebot/old
newpage=~/.firebot/newpage
olddata=~/.firebot/old_data
newdata=~/.firebot/fds_times.csv
running=~/.fdssmvgit/bot_running
curdir=`pwd`
EXIT="yes"

# don't update status page if firebot is running
cd $firebotdir
if [ -e $running ] ; then
  exit
fi
./status2html.sh  > $new
if [ -e $old ]; then
  ndiff=`diff $old $new|wc -l`
  if [ "$ndiff" == "0" ] ; then
    EXIT="no"
  else
    cp $new $old
  fi
fi
if [ -e $olddata ]; then
  ndiff=`diff $olddata $newdata|wc -l`
  if [ "$ndiff" == "0" ] ; then
     EXIT="no"
  else
     cp $newdata $olddata
  fi
fi
if [ "$EXIT" == "yes" ]; then
   exit
fi

./status_pubtop.sh > $newpage
./status2html.sh >> $newpage
./status_pubbot.sh >> $newpage

cd $gitwebrepo
git remote update
git merge origin/gh-pages

cp $newpage firebot_status.html
git add firebot_status.html

cp $newdata fds_times.csv
git add fds_times.csv

git commit -m "firebot: update firebot status page `date`"
git push
