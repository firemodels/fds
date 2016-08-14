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
RUNSMOKEBOT=1
MOVIE=
SSH=
UPLOAD=
FORCE=
COMPILER=intel
SMOKEBOT_LITE=
TESTFLAG=

WEB_URL=
web_DIR=/var/www/html/`whoami`
if [ -d $web_DIR ]; then
  IP=`wget http://ipinfo.io/ip -qO -`
  HOST=`host $IP | awk '{printf("%s\n",$5);}'`
  WEB_URL=http://$HOST/`whoami`
else
  web_DIR=
fi

# checking to see if a queing system is available
QUEUE=smokebot
notfound=`qstat -a 2>&1 | tail -1 | grep "not found" | wc -l`
if [ $notfound -eq 1 ] ; then
  QUEUE=none
fi

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
echo "-I compiler - intel or gnu [default: $COMPILER]"
if [ "$EMAIL" != "" ]; then
echo "-m email_address - [default: $EMAIL]"
else
echo "-m email_address"
fi
echo "-q queue [default: $QUEUE]"
echo "-L - smokebot lite,  run only stages that build a debug fds and run cases with it"
echo "                    (no release fds, no release cases, no manuals, etc)"
echo "-M  - make movies"
echo "-r - FDS-SMV repository location [default: $FDSREPO]"
echo "-S host - generate images on host"
echo "-t - use test smokeview"
echo "-u - update repo"
echo "-U - upload guides"
echo "-v - show options used to run smokebot"
if [ "$web_DIR" == "" ]; then
echo "-w directory - web directory containing summary pages"
else
echo "-w directory - web directory containing summary pages [default: $web_DIR]"
fi
if [ "$WEB_URL" == "" ]; then
echo "-W url - web url of summary pages"
else
echo "-W url - web url of summary pages [default: $WEB_URL]"
fi
exit
}

while getopts 'aAb:C:cd:fhI:Lm:Mq:r:S:tuUvw:W:' OPTION
do
case $OPTION  in
  a)
   RUNAUTO=-a
   ;;
  A)
   RUNAUTO=-A
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
  L)
   SMOKEBOT_LITE="-L"
   ;;
  m)
   EMAIL="$OPTARG"
   ;;
  M)
   MOVIE="-M"
   ;;
  q)
   QUEUE="$OPTARG"
   ;;
  r)
   FDSREPO="$OPTARG"
   ;;
  S)
   SSH="-S $OPTARG"
   ;;
  t)
   TESTFLAG="-t"
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
  w)
   web_DIR="$OPTARG"
   ;;
  W)
   WEB_URL="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

if [ ! "$web_DIR" == "" ]; then
  web_DIR="-w $web_DIR"
fi
if [ ! "$WEB_URL" == "" ]; then
  WEB_URL="-W $WEB_URL"
fi

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

QUEUE="-q $QUEUE"

if [ "$EMAIL" != "" ]; then
  EMAIL="-m $EMAIL"
fi

if [[ "$RUNSMOKEBOT" == "1" ]]; then
  if [[ "$UPDATEREPO" == "-u" ]]; then
     cd $FDSREPO
     git remote update &> /dev/null
     git checkout $BRANCH &> /dev/null
     git merge origin/$BRANCH &> /dev/null
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
  ./$botscript $TESTFLAG $RUNAUTO $COMPILER $SSH $SMOKEBOT_LITE $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $web_DIR $WEB_URL $UPDATEREPO $QUEUE $UPLOAD $EMAIL $MOVIE "$@"
  rm $running
else
  echo ./$botscript $TESTFLAG $RUNAUTO $COMPILER $SMOKEBOT_LITE $SSH $BRANCH $CFASTREPO $FDSREPO $CLEANREPO $web_DIR $WEB_URL $UPDATEREPO $QUEUE $UPLOAD $EMAIL $MOVIE "$@"
fi
