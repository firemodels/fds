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
cd $CURDIR/cor
git remote add firemodels git@github.com:firemodels/cor.git
git remote update

cd $CURDIR/exp
git remote add firemodels git@github.com:firemodels/exp.git
git remote update

cd $CURDIR/out
git remote add firemodels git@github.com:firemodels/out.git
git remote update

cd $CURDIR/smv
git remote add firemodels git@github.com:firemodels/smv.git
git remote update

cd $CURDIR/fds
git remote add firemodels git@github.com:firemodels/fds.git
git remote update

cd $CURDIR/radcal
git remote add firemodels git@github.com:firemodels/radcal.git
git remote update

cd $CURDIR
