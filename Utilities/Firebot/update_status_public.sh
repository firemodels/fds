#!/bin/bash
repo=~/FDS-SMVgitclean
firebotdir=~/firebotgit
old=$firebotdir/history/old
new=$firebotdir/history/new
running=$firebotdir/bot_running
curdir=`pwd`

# don't update status page if firebot is running
cd $firebotdir
if [ -e $running ] ; then
  exit
fi
./list2html.sh statusonly > $new
ndiff=`diff $old $new|wc -l`
if [ "$ndiff" == "0" ] ; then
exit
fi
cp $new $old
cd $repo
git remote update
git checkout gh-pages
cp $new firebot_status.html
git add firebot_status.html
git commit -m "firebot: update firebot status page `date`"
git push

