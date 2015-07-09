#!/bin/bash

CURDIR=`pwd`
gitrepo=FDS-SMVgitclean
gitrepodir=~/$gitrepo
botdir=~/firebotgit

if [ ! -d $gitrepodir ] ; then
  cd 
  echo $gitrepodir does not exist - creating
  git clone git@github.com:firemodels/fds-smv.git $gitrepo
  echo $gitrepodir created.
fi

# create directory where firebot runs

if [ ! -d $botdir ] ; then
  cd 
  echo $botdir does not exist - creating
  mkdir $botdir
  cd %botdir%
  cp $gitrepodir/Utilities/Firebot/*.sh .
  echo $botdir created.
fi

gitrepo=cfastgitclean
gitrepodir=~/$gitrepo
if [ ! -d $gitrepodir ] ; then
  cd 
  echo $gitrepodir does not exist - creating
  git clone git@github.com:firemodels/cfast.git $gitrepo
  echo $gitrepodir created.
fi

cd $CURDIR
