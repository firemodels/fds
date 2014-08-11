#!/bin/bash  

MAKEMOVIES=
while getopts 'm' OPTION
do case $OPTION in
  m)
   MAKEMOVIES="-m"
  ;;
esac
done
bash -lc "./smokebot_linux_wrapper.sh -q smokebot $MAKEMOVIES" &
