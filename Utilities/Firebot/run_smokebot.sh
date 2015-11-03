#!/bin/bash
running=smokebot_running

CURDIR=`pwd`
FDSREPO=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  FDSREPO=$FDSSMV
fi
if [ -e .fds_git ]; then
  cd ../..
  FDSREPO=`pwd`
  cd $CURDIR
fi
CFASTREPO=~/cfastgitclean

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
UPLOAD=

function usage {
echo "run_smokebot.sh -a -b branch_name -c -h -m email_address -r repo location -S host -u -v"
echo "Run smokebot_linux.sh V&V testing script"
echo ""
echo "Options:"
echo "-a - run automatically if FDS or smokeview source has changed"
echo "-b - branch_name - run smokebot using the branch branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-C - cfast repository location [default: $CFASTREPO]"
echo "-h - display this message"
echo "-m email_address"
echo "-q queue"
echo "-M  - make movies"
echo "-r - FDS-SMV repository location [default: $FDSREPO]"
echo "-S host - generate images on host"
echo "-u - update repo"
echo "-U - upload guides"
echo "-v - show options used to run smokebot"
exit
}

while getopts 'ab:C:cd:hm:Mq:r:S:uUv' OPTION
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
  C)
   CFASTREPO="-C $OPTARG"
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
   FDSREPO="$OPTARG"
   ;;
  S)
   SSH="-S $OPTARG"
   ;;
  u)
   UPDATEREPO=-u
   ;;
  U)
   UPLOAD="-U"
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

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [[ "$UPDATEREPO" == "-u" ]]; then
     cd $FDSREPO
     git remote update
     git checkout $BRANCH
     git pull
     cp Utilities/Firebot/$botscript $CURDIR/.
     cd $CURDIR
  fi
fi
CFASTREPO="-C $CFASTREPO"
FDSREPO="-r $FDSREPO"
BRANCH="-b $BRANCH"
if [[ "$RUNSMOKEBOT" == "1" ]]; then
  touch $running
  ./$botscript $RUNAUTO $SSH $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $UPDATEREPO $QUEUE $UPLOAD $MAILTO $MOVIE "$@"
  rm $running
else
  echo ./$botscript $RUNAUTO $SSH $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $UPDATEREPO $QUEUE $UPLOAD $MAILTO $MOVIE "$@"
fi
