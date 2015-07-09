#!/bin/bash
if [ -e bot_running ] ; then
  exit
fi

CURDIR=`pwd`
FDS_GITbase=FDS-SMVgitclean
BRANCH=development
botscript=firebot_linux.sh
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

if [ "$FDS_GITbase" != "" ] ; then
   cFDS_GITbase=-d $FDS_GITbase
fi 
if [ "$BRANCH" != "" ] ; then
   cFDS_GITbase=-b $BRANCH
fi 
if [ "$UPDATEREPO" == "1" ] ; then
   cd ~/$FDS_GITbase
   git checkout $BRANCH
   git pull
   cp Utilities/Firebot/$botscript $CURDIR/.
   cd $CURDIR
fi
touch bot_running
./$botscript $cBRANCH $cFDS_GITbase "$@"
rm bot_running
