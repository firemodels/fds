#!/bin/bash
CURDIR=`pwd`
repo=FDS-SMVgitclean
while getopts 'r:' OPTION
do
case $OPTION  in
  r)
   repo="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))


cd ~/$repo
git remote update
git pull
cd $CURDIR
cp ~/$repo/Utilities/Firebot/*.sh .
