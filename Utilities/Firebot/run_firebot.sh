#!/bin/bash

running=bot_running
if [ -e bot_running ] ; then
  echo Firebot is already running.
  echo Erase the file $running if this is not the case.
  exit
fi

CURDIR=`pwd`
FDS_GITbase=
FDS_PARENT=
BRANCH=development
botscript=firebot_linux.sh
cFDS_GITbase=
cBRANCH=
UPDATEREPO=
while getopts 'b:d:p:u' OPTION
do
case $OPTION  in
  b)
   BRANCH="$OPTARG"
   ;;
  d)
   FDS_GITbase="$OPTARG"
   ;;
  p)
   FDS_PARENT="$OPTARG"
   ;;
  u)
   UPDATEREPO=1
   ;;
esac
done
shift $(($OPTIND-1))

if [[ "$FDS_GITbase" != "" ]]; then
   FDS_GITbase="-d $FDS_GITbase"
fi 
if [[ "$FDS_PARENT" != "" ]]; then
   FDS_PARENT="-p $FDS_PARENT"
fi 
if [[ "$BRANCH" != "" ]]; then
   cBRANCH="-b $BRANCH"
fi 
if [[ "$UPDATEREPO" == "1" ]]; then
   cd ~/$FDS_GITbase
   git checkout $BRANCH
   git pull
   cp Utilities/Firebot/$botscript $CURDIR/.
   cd $CURDIR
fi
touch $running
./$botscript $cBRANCH $FDS_GITbase $FDS_PARENT "$@"
rm $running
