#!/bin/bash
running=bot_running
if [ -e bot_running ] ; then
  echo Firebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi

reponame=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  reponame=$FDSSMV
fi

function usage {
echo "run_firebot.sh [ -b branch_name -h -r repo_name -u ]"
echo "Run Firebot V&V testing script"
echo ""
echo "Options:"
echo "-b - branch_name - run firebot using branch_name [default: $BRANCH]"
echo "-h - display this message"
echo "-r - repository location [default: $reponame]"
echo "-u - update repo"
echo "-v - show options used to run firebot"
exit
}

CURDIR=`pwd`
BRANCH=development
botscript=firebot_linux.sh
UPDATEREPO=
RUNFIREBOT=1
while getopts 'b:hr:uv' OPTION
do
case $OPTION  in
  b)
   BRANCH="$OPTARG"
   ;;
  h)
   usage;
   ;;
  r)
   reponame="$OPTARG"
   ;;
  u)
   UPDATEREPO=1
   ;;
  v)
   RUNFIREBOT=
   ;;
esac
done
shift $(($OPTIND-1))
exit
if [[ "$UPDATEREPO" == "1" ]]; then
   cd $reponame
   git remote update
   git checkout $BRANCH
   git pull
   cp Utilities/Firebot/$botscript $CURDIR/.
   cd $CURDIR
fi
touch $running
BRANCH="-b $BRANCH"
reponame="-r $reponame"
if [ "$RUNFIREBOT" == "1" ] ; then
  ./$botscript $BRANCH $reponame "$@"
else
  echo ./$botscript $BRANCH $reponame "$@"
fi
rm $running
