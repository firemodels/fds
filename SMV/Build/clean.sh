#!/bin/bash
clean=1
while getopts 'n' OPTION
do
case $OPTION in
  n)
   clean=0
   ;;
esac
done
if [ "$clean" == "1" ]; then
  echo cleaning files
  rm -f *.o
else
  echo not cleaning files
fi
