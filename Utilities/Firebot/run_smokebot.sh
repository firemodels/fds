#!/bin/bash

running=bot_running
if [ -e bot_running ] ; then
  echo Smokebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi

CURDIR=`pwd`
FDS_GITbase=FDS-SMVgitclean
BRANCH=development
botscript=smokebot_linux.sh
cFDS_GITbase=
cBRANCH=
UPDATEREPO=
while getopts 'b:d:u' OPTION
do
case $OPTION  in
  b)
   BRANCH="$OPTARG"
   ;;
  d)
   FDS_GITbase="$OPTARG"
   ;;
  u)
   UPDATEREPO=1
   ;;
esac
done
shift $(($OPTIND-1))

if [[ "$FDS_GITbase" != "" ]]; then
   cFDS_GITbase="-d $FDS_GITbase"
fi 
if [[ "$BRANCH" != "" ]]; then
   cFDS_GITbase="-b $BRANCH"
fi 
if [[ "$UPDATEREPO" == "1" ]]; then
   cd ~/$FDS_GITbase
   git checkout $BRANCH
   git pull
   cp Utilities/Firebot/$botscript $CURDIR/.
   cd $CURDIR
fi
touch $running
./$botscript $cBRANCH $cFDS_GITbase "$@"
rm $running
