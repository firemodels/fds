#!/bin/bash

FDSREPO=~/FDS-SMVgitclean
if [ "$FDSSMV" != "" ] ; then
  FDSREPO=$FDSSMV
fi
function usage {
echo "update submodules in a git repo. This script assumes"
echo "that the command"
echo ""
echo "git submodule update --init --recursive"
echo ""
echo "was run initially to set up the submodule."
echo ""
echo "Options:"
echo "-h - display this message"
echo "-r - repository location [default: $FDSREPO]"
exit
}
while getopts 'hr:' OPTION
do
case $OPTION  in
  h)
   usage
   exit
   ;;
  r)
   FDSREPO="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

curdir=`pwd`
cd $FDSREPO
echo updating submodules in $FDSREPO
git submodule foreach git remote update 
git submodule foreach git merge origin/master
cd $curdir
