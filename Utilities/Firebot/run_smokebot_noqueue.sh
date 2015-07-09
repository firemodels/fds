#!/bin/bash  

MAKEMOVIES=
while getopts 'm' OPTION
do case $OPTION in
  m)
   MAKEMOVIES="-m"
  ;;
esac
done
bash -lc "./run_smokebot.sh -s -q none $MAKEMOVIES" &
