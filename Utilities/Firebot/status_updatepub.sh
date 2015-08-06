#!/bin/bash
repo=~/FDS-SMVgitclean
firebotdir=~/firebotgit
old=$firebotdir/history/old
new=$firebotdir/history/new
newpage=$firebotdir/history/newpage
running=$firebotdir/bot_running
curdir=`pwd`

# don't update status page if firebot is running
cd $firebotdir
if [ -e $running ] ; then
  exit
fi
./status2html.sh  > $new
ndiff=`diff $old $new|wc -l`
if [ "$ndiff" == "0" ] ; then
exit
fi
./status_pubtop.sh > $newpage
./status2html.sh >> $newpage
./status_pubbot.sh >> $newpage
cp $new $old
cd $repo
git remote update
git checkout gh-pages
cp $newpage firebot_status.html
git add firebot_status.html
git commit -m "firebot: update firebot status page `date`"
git push
git checkout development
