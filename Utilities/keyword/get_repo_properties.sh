#!/bin/bash
repo_dir=$1

export revision=unknown
export revision_date=unknown
export build_date=unknown
validsvn=0
validgit=0
havesvn=1
havegit=1

CURDIR=`pwd`

# checking if repo_dir exists

if [ ! -e "$repo_dir" ] ; then
  echo directory $repo_dir does not exist
  exit
fi

# looking for svn

notfound=`svn info 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ] ; then
  havesvn=0
fi

# looking for git

notfound=`git 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ] ; then
  havegit=0
fi

# warn if both svn and git are not available

if [ "$havesvn" == "0" ] ; then
  if [ "$havesvn" == "0" ] ; then
    echo "*** warning: neither svn nor git are available"
    cd $CURDIR
    exit
  fi
fi

# check if repo_dir is a valid svn repository

if [ "$havesvn" == "1" ]; then
  validsvn=1
  notworking=`svn info 2>&1 | grep "not a working copy" | wc -l`
  if [ "$notworking" == "1" ] ; then
    validsvn=0
  fi
fi

# if not an svn repo check if repo_dir is a valid git repository

if [ "$validsvn" == "0" ] ; then
  if [ "$havegit" == "1" ] ; then
    validgit=1
    notworking=`git log . 2>&1 | grep "Not a git repository" | wc -l`
    if [ "$notworking" == "1" ] ; then
      validgit=0
    fi
  fi
fi

# warn if not a valid repository

if [ "$validsvn" == "0" ] ; then
  if [ "$validgit" == "0" ] ; then
    echo "*** warning: $repo_dir is not a valid repository (git or svn)"
    cd $CURDIR
    exit
  fi
fi

# examine properties of repo_dir 

cd $repo_dir

# get revision number

if [ "$validsvn" == 1 ] ; then
  revision=`svn info 2>&1 | grep "Last Changed Rev:" | awk -F' ' '{print $4}'`
fi
if [ "$validgit" == 1 ] ; then
  revision=`git log . 2>&1 | head -1 | awk -F " " '{print $2}'`
fi

# get date/time of latest repository commit

if [ "$validsvn" == 1 ] ; then
  revision_date=`svn info 2>&1 | grep "Last Changed Date:" | awk -F" " '{print $4,$5}'`
fi
if [ "$validgit" == 1 ] ; then
  revision_date=`git log --date=iso . 2>&1 | head -3 | tail -1 | awk -F" " '{print $2,$3}'`
fi

# get current date/time

build_date=`date "+%F %T"`

cd $CURDIR
export revision
export revision_date
export build_date
