#!/bin/bash
gitwebrepo=~/FDS-SMVgitweb
firebotdir=~/FDS-SMVgitclean/Utilities/Firebot
old=~/.firebot/old
new=~/.firebot/new
newpage=~/.firebot/newpage
running=~/.fdssmvgit/bot_running
curdir=`pwd`

# don't update status page if firebot is running
cd $firebotdir
if [ -e $running ] ; then
  exit
fi
./status2html.sh  > $new
if [ -e $old ]; then
  ndiff=`diff $old $new|wc -l`
  if [ "$ndiff" == "0" ] ; then
    exit
  fi
fi
./status_pubtop.sh > $newpage
./status2html.sh >> $newpage
./status_pubbot.sh >> $newpage
cp $new $old
cd $gitwebrepo
git remote update
git merge origin/gh-pages
cp $newpage firebot_status.html
git add firebot_status.html
git commit -m "firebot: update firebot status page `date`"
git push
