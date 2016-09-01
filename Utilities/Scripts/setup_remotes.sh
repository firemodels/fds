#!/bin/bash

function usage {
echo "Setup remote connections with each firemodels repo"
echo ""
echo "Options:"
echo "-h - display this message"
exit
}

while getopts 'h' OPTION
do
case $OPTION  in
  h)
   usage;
   ;;
esac
done
shift $(($OPTIND-1))

echo "You are about to setup remotes for the repos:"
echo "cor, exp, out, smv, fds and radcal at "
echo "git@github.com:firemodels"
echo "Press any key to continue or <CTRL> c to abort."
echo "You should only run this command if you have cloned these" 
echo "repos from repos that have been forked from firemodels"
read val

CURDIR=`pwd`
for repo in cfast cor exp fds out radcal smv
do
echo.
echo updating remotes for $repo
  cd $CURDIR/$repo
  git remote add firemodels git@github.com:firemodels/$repo.git
  git remote update
done

