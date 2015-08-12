#!/bin/bash

running=bot_running

reponame=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  reponame=$FDSSMV
fi

CURDIR=`pwd`
FDS_GITbase=FDS-SMVgitclean
BRANCH=development
botscript=smokebot_linux.sh
RUNAUTO=
UPDATEREPO=
QUEUE=
RUNSMOKEBOT=1
fopt=
mopt=
while getopts 'ab:d:fmq:r:uv' OPTION
do
case $OPTION  in
  a)
   RUNAUTO=-a
   ;;
  b)
   BRANCH="$OPTARG"
   ;;
  f)
   fopt="-f"
   ;;
  m)
   mopt="-m"
   ;;
  r)
   reponame="$OPTARG"
   ;;
  u)
   UPDATEREPO=1
   ;;
  v)
   RUNSMOKEBOT=
   ;;
esac
done
shift $(($OPTIND-1))

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [ -e bot_running ] ; then
    echo Smokebot is already running.
    echo Erase the file $running if this is not the case.
    exit
  fi
fi

FDS_GITBASE=`basename $reponame`
if [[ "$QUEUE" != "" ]]; then
   QUEUE="-q $QUEUE"
fi 
reponame="-r $reponame"
if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [[ "$UPDATEREPO" == "1" ]]; then
     cd ~/$FDS_GITBASE
     git remote update
     git checkout $BRANCH
     git pull
     cp Utilities/Firebot/$botscript $CURDIR/.
     cd $CURDIR
  fi
fi
BRANCH="-b $BRANCH"
if [[ "$RUNSMOKEBOT" == "1" ]]; then
  touch $running
  ./$botscript $RUNAUTO $BRANCH $reponame $QUEUE $fopt $mopt "$@"
  rm $running
else
  echo ./$botscript $RUNAUTO $BRANCH $reponame $QUEUE $fopt $mopt "$@"
fi
