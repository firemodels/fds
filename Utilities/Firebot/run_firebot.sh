#!/bin/bash
if [ ! -d ~/.fdssmvgit ] ; then
  mkdir ~/.fdssmvgit
fi
firebot_pid=~/.fdssmvgit/firebot_pid

CURDIR=`pwd`
reponame=~/FDS-SMVgitclean
if [ "$FIREMODELS" != "" ] ; then
  reponame=$FIREMODELS
fi
if [ -e .fds_git ]; then
  cd ../../..
  reponame=`pwd`
  cd $CURDIR
else
  echo "***error: firebot not running in the FDS repo"
  exit
fi

# checking to see if a queing system is available
QUEUE=firebot
notfound=`qstat -a 2>&1 | tail -1 | grep "not found" | wc -l`
if [ $notfound -eq 1 ] ; then
  QUEUE=none
fi

function usage {
echo "Verification and validation testing script for FDS"
echo ""
echo "Options:"
echo "-b - branch_name - run firebot using branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-f - force firebot run"
echo "-F - skip figure generation and build document stages"
echo "-h - display this message"
echo "-i - use installed version of smokeview"
echo "-k - kill firebot if it is running"
echo "-L - firebot lite,  run only stages that build a debug fds and run cases with it"
echo "                    (no release fds, no release cases, no matlab, etc)"
if [ "$EMAIL" != "" ]; then
echo "-m email_address [default: $EMAIL]"
else
echo "-m email_address "
fi
echo "-q queue - specify queue [default: $QUEUE]"
echo "-r - repository location [default: $reponame]"
echo "-s - skip matlab and build document stages"
echo "-u - update repo"
echo "-U - upload guides (only by user firebot)"
echo "-v - show options used to run firebot"
exit
}

LIST_DESCENDANTS ()
{
  local children=$(ps -o pid= --ppid "$1")

  for pid in $children
  do
    LIST_DESCENDANTS "$pid"
  done

  echo "$children"
}

USEINSTALL=
BRANCH=master
botscript=firebot.sh
UPDATEREPO=
CLEANREPO=0
UPDATE=
CLEAN=
RUNFIREBOT=1
UPLOADGUIDES=
FORCE=
SKIPMATLAB=
SKIPFIGURES=
FIREBOT_LITE=
KILL_FIREBOT=
while getopts 'b:cFfhikLm:q:nr:suUv' OPTION
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
  F)
   SKIPFIGURES=-F
   ;;
  h)
   usage;
   ;;
  i)
   USEINSTALL="-i"
   ;;
  k)
   KILL_FIREBOT="1"
   ;;
  L)
   FIREBOT_LITE=-L
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
  s)
   SKIPMATLAB=-s
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

if [ "$KILL_FIREBOT" == "1" ]; then
  if [ -e $firebot_pid ] ; then
    PID=`head -1 $firebot_pid`
    echo killing process invoked by firebot
    kill -9 $(LIST_DESCENDANTS $PID)
    echo "killing firebot (PID=$PID)"
    kill -9 $PID
    JOBIDS=`qstat -a | grep FB_ | awk -v user="$USER" '{if($2==user){print $1}}'`
    if [ "$JOBIDS" != "" ]; then
      echo killing firebot jobs with Id:$JOBIDS
      qdel $JOBIDS
    fi
    echo firebot process $PID killed
    if [ -e $firebot_pid ]; then
      rm $firebot_pid
    fi
  else
    echo firebot is not running, cannot be killed.
  fi
  exit
fi
if [ -e $firebot_pid ] ; then
  if [ "$FORCE" == "" ] ; then
    echo Firebot or smokebot are already running. If this
    echo "is not the case re-run using the -f option."
    exit
  fi
fi
if [[ "$EMAIL" != "" ]]; then
  EMAIL="-m $EMAIL"
fi
if [[ "$UPDATEREPO" == "1" ]]; then
   UPDATE=-u
   cd $reponame/fds
   if [[ "$RUNFIREBOT" == "1" ]]; then
     git fetch origin &> /dev/null
     git checkout $BRANCH &> /dev/null
     git merge origin/$BRANCH &> /dev/null
     cd Utilities/Firebot
     FIREBOTDIR=`pwd`
     if [[ "$CURDIR" != "$FIREBOTDIR" ]]; then
        echo "***error: firebot not running in the $FIREBOTDIR"
        exit
     fi
     cd $CURDIR
  fi
fi
if [[ "$CLEANREPO" == "1" ]]; then
  CLEAN=-c
fi
BRANCH="-b $BRANCH"
QUEUE="-q $QUEUE"
reponame="-r $reponame"
if [ "$RUNFIREBOT" == "1" ] ; then
  touch $firebot_pid
  ./$botscript -p $firebot_pid $UPDATE $FIREBOT_LITE $USEINSTALL $UPLOADGUIDES $CLEAN $BRANCH $QUEUE $SKIPMATLAB $SKIPFIGURES $reponame $EMAIL "$@"
else
  echo ./$botscript $FIREBOT_LITE $UPDATE $USEINSTALL $UPLOADGUIDES $CLEAN $BRANCH $QUEUE $SKIPMATLAB $SKIPFIGURES $reponame $EMAIL "$@"
fi
rm $firebot_pid
