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
MOVIE=
SSH=
MAILTO=

function usage {
echo "run_smokebot.sh -a -b branch_name -c -h -m email_address -r repo location -S host -u -v"
echo "Run smokebot_linux.sh V&V testing script"
echo ""
echo "Options:"
echo "-a - run automatically if FDS or smokeview source has changed"
echo "-b - branch_name - run smokebot using the branch branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-h - display this message"
echo "-m email_address"
echo "-q queue"
echo "-M  - make movies"
echo "-r - repository location [default: $reponame]"
echo "-S host - generate images on host"
echo "-u - update repo"
echo "-v - show options used to run smokebot"
exit
}

while getopts 'ab:cd:hm:Mq:r:S:uv' OPTION
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
  h)
   usage
   exit
   ;;
  m)
   MAILTO="-m $OPTARG"
   ;;
  M)
   MOVIE="-M"
   ;;
  q)
   QUEUE="-q $OPTARG"
   ;;
  r)
   reponame="$OPTARG"
   ;;
  S)
   SSH="-S $OPTARG"
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
  ./$botscript $RUNAUTO $SSH $BRANCH $reponame $CLEANREPO $UPDATEREPO $QUEUE $MAILTO $MOVIE "$@"
  rm $running
else
  echo ./$botscript $RUNAUTO $SSH $BRANCH $reponame $CLEANREPO $UPDATEREPO $QUEUE $MAILTO $MOVIE "$@"
fi
