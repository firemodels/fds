#!/bin/bash
UN=$USER

function usage {
echo "Download firemodels repos from github"
echo ""
echo "Options:"
echo "-u user - specify user [default: $UN]"
echo "          specify firemodels to download central repos"
echo "-h - display this message"
exit
}

while getopts 'hu:' OPTION
do
case $OPTION  in
  h)
   usage;
   ;;
  u)
   UN="$OPTARG"
   ;;
esac
done
shift $(($OPTIND-1))

echo "You are about to download the repos:"
echo "cor, exp, out, smv, fds and radcal from "
echo "git@github.com:$UN"
echo "Press any key to continue, <CTRL> c to abort or"
echo "$0 -h"
echo "for other options"
read val

for repo in cor cfast exp out radcal smv
do 
  echo
  git clone git@github.com\:$UN/$repo.git
done
echo
git clone  --recursive git@github.com\:$UN/fds.git
