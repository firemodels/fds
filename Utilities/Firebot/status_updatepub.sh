#!/bin/bash
repo=~/FDS-SMVgitweb
firebotdir=~/firebotgit
old=$firebotdir/history/old
new=$firebotdir/history/new
newpage=$firebotdir/history/newpage
running=$firebotdir/firebot_running
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
git pull
cp $newpage firebot_status.html
git add firebot_status.html
git commit -m "firebot: update firebot status page `date`"
git push
