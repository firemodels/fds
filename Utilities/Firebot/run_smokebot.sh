#!/bin/bash
if [ ! -d ~/.fdssmvgit ] ; then
  mkdir ~/.fdssmvgit
fi
running=~/.fdssmvgit/bot_running

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
botscript=smokebot.sh
RUNAUTO=
CLEANREPO=
UPDATEREPO=
QUEUE=
RUNSMOKEBOT=1
MOVIE=
SSH=
MAILTO=
UPLOAD=
FORCE=
COMPILER=intel

function usage {
echo "Verification and validation testing script for smokeview"
echo ""
echo "Options:"
echo "-a - run automatically if FDS or smokeview source has changed"
echo "-b - branch_name - run smokebot using the branch branch_name [default: $BRANCH]"
echo "-c - clean repo"
echo "-C - cfast repository location [default: $CFASTREPO]"
echo "-f - force smokebot run"
echo "-h - display this message"
echo "-I - specify compiler (intel or gnu)"
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

while getopts 'ab:C:cd:fhI:m:Mq:r:S:uUv' OPTION
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
   CFASTREPO="$OPTARG"
   ;;
  I)
   COMPILER="$OPTARG"
   ;;
  f)
   FORCE=1
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

COMPILER="-I $COMPILER"

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [ "$FORCE" == "" ]; then
    if [ -e $running ] ; then
      echo Smokebot or firebot are already running.
      echo "Re-run using the -f option if this is not the case."
      exit
    fi
  fi
fi

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [[ "$UPDATEREPO" == "-u" ]]; then
     cd $FDSREPO
     git remote update &> /dev/null
     git checkout $BRANCH &> /dev/null
     git pull &> /dev/null
     cd Utilities/Firebot
     FIREBOTDIR=`pwd`
     if [ "$FIREBOTDIR" != "$CURDIR" ]; then
       cp $botscript $CURDIR/.
     fi
     cd $CURDIR
  fi
fi
CFASTREPO="-C $CFASTREPO"
FDSREPO="-r $FDSREPO"
BRANCH="-b $BRANCH"
if [[ "$RUNSMOKEBOT" == "1" ]]; then
  touch $running
  ./$botscript $RUNAUTO $COMPILER $SSH $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $UPDATEREPO $QUEUE $UPLOAD $MAILTO $MOVIE "$@"
  rm $running
else
  echo ./$botscript $RUNAUTO $COMPILER $SSH $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $UPDATEREPO $QUEUE $UPLOAD $MAILTO $MOVIE "$@"
fi
