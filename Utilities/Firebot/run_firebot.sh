#!/bin/bash
running=firebot_running
if [ -e $running ] ; then
  echo Firebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi

QUEUE=firebot
reponame=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  reponame=$FDSSMV
fi

function usage {
echo "run_firebot.sh [ -b branch_name -c -h -m email_address -r repo location -u -v]"
echo "Run Firebot V&V testing script"
echo ""
echo "Options:"
echo "-b - branch_name - run firebot using branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-h - display this message"
echo "-m email_address "
echo "-q - queue_name - run cases using the queue queue_name"
echo "     default: $QUEUE"
echo ""
echo "-r - repository location [default: $reponame]"
echo "-u - update repo"
echo "-U - upload guides (only by user firebot)"
echo "-v - show options used to run firebot"
exit
}

CURDIR=`pwd`
BRANCH=development
botscript=firebot_linux.sh
UPDATEREPO=
CLEANREPO=0
UPDATE=
CLEAN=
RUNFIREBOT=1
EMAIL=
UPLOADGUIDES=
while getopts 'b:chm:q:nr:uUv' OPTION
do
case $OPTION  in
  b)
   BRANCH="$OPTARG"
   ;;
  c)
   CLEANREPO=1
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
  ./$botscript $UPDATE $UPLOADGUIDES $CLEAN $BRANCH $QUEUE $reponame $EMAIL "$@"
else
  echo ./$botscript $UPDATE $UPLOADGUIDES $CLEAN $BRANCH $QUEUE $reponame $EMAIL "$@"
fi
rm $running
