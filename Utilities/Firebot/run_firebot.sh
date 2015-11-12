#!/bin/bash
running=firebot_running

CURDIR=`pwd`
QUEUE=firebot
reponame=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  reponame=$FDSSMV
fi
if [ -e .fds_git ]; then
  cd ../..
  reponame=`pwd`
  cd $CURDIR
fi


function usage {
echo "Verification and validation testing script for FDS"
echo ""
echo "Options:"
echo "-b - branch_name - run firebot using branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-f - force firebot run"
echo "-h - display this message"
echo "-m email_address "
echo "-q - queue_name - run cases using the queue queue_name"
echo "     default: $QUEUE"
echo "-r - repository location [default: $reponame]"
echo "-S host - generate images on host"
echo "-u - update repo"
echo "-U - upload guides (only by user firebot)"
echo "-v - show options used to run firebot"
exit
}

BRANCH=development
botscript=firebot_linux.sh
UPDATEREPO=
CLEANREPO=0
UPDATE=
CLEAN=
RUNFIREBOT=1
EMAIL=
UPLOADGUIDES=
SSH=
FORCE=
while getopts 'b:cfhm:q:nr:S:uUv' OPTION
do
case $OPTION  in
  b)
   BRANCH="$OPTARG"
   ;;
  c)
   CLEANREPO=1
   ;;
  f)
   FORCE=1
   ;;
  h)
   usage;
   ;;
  m)
   EMAIL="$OPTARG"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  n)
   UPDATEREPO=0
   ;;
  r)
   reponame="$OPTARG"
   ;;
  S)
   SSH="-S $OPTARG"
   ;;
  u)
   UPDATEREPO=1
   ;;
  U)
   UPLOADGUIDES=-U
   ;;
  v)
   RUNFIREBOT=0
   ;;
esac
done
shift $(($OPTIND-1))

if [ -e $running ] ; then
  if [ "$FORCE" == ""] ; then
    echo Firebot is already running.
    echo Erase the file $running if this is not the case
    echo or rerun using the -f option.
    exit
  fi
fi
if [[ "$EMAIL" != "" ]]; then
  EMAIL="-m $EMAIL"
fi
if [[ "$UPDATEREPO" == "1" ]]; then
   UPDATE=-u
   cd $reponame
   if [[ "$RUNFIREBOT" == "1" ]]; then
     git remote update
     git checkout $BRANCH
     git merge origin/$BRANCH
     cd Utilities/Firebot
     FIREBOTDIR=`pwd`
     if [[ "$CURDIR" != "$FIREBOTDIR" ]]; then
       cp $botscript $CURDIR/.
     fi
     cd $CURDIR
  fi
fi
if [[ "$CLEANREPO" == "1" ]]; then
  CLEAN=-c
fi
touch $running
BRANCH="-b $BRANCH"
QUEUE="-q $QUEUE"
reponame="-r $reponame"
if [ "$RUNFIREBOT" == "1" ] ; then
  ./$botscript $UPDATE $UPLOADGUIDES $SSH $CLEAN $BRANCH $QUEUE $reponame $EMAIL "$@"
else
  echo ./$botscript $UPDATE $UPLOADGUIDES $SSH $CLEAN $BRANCH $QUEUE $reponame $EMAIL "$@"
fi
rm $running
