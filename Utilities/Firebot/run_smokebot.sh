#!/bin/bash
running=smokebot_running

reponame=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  reponame=$FDSSMV
fi

CURDIR=`pwd`
FDS_GITbase=FDS-SMVgitclean
BRANCH=development
botscript=smokebot_linux.sh
RUNAUTO=
CLEANREPO=
UPDATEREPO=
QUEUE=
RUNSMOKEBOT=1
fopt=
mopt=
while getopts 'ab:cd:fmq:r:uv' OPTION
do
case $OPTION  in
  a)
   RUNAUTO=-a
   ;;
  b)
   BRANCH="$OPTARG"
   ;;
  c)
   CLEANREPO=-c
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
   UPDATEREPO=-u
   ;;
  v)
   RUNSMOKEBOT=
   ;;
esac
done
shift $(($OPTIND-1))

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [ -e $running ] ; then
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
  if [[ "$UPDATEREPO" == "-u" ]]; then
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
  ./$botscript $RUNAUTO $BRANCH $reponame $CLEANREPO $UPDATEREPO $QUEUE $fopt $mopt "$@"
  rm $running
else
  echo ./$botscript $RUNAUTO $BRANCH $reponame $CLEANREPO $UPDATEREPO $QUEUE $fopt $mopt "$@"
fi
