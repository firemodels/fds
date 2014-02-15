#!/bin/bash
SIZE="-m64"
COMPILER="icc"
COMPILER2="icc"
PLATFORM=""
while getopts '36ghit' OPTION
do
case $OPTION in
  h)
  echo "options:"
  echo "-3 - build a 32 bit version"
  echo "-t - build a 32 bit version but don't use -m32 option"
  echo "-6 - build a 64 bit version"
  echo "-g - use the gcc compiler"
  echo "-i - use the Intel icc compiler"
  exit
  ;;
  3)
   SIZE="-m32"
  ;;
  6)
   SIZE="-m64"
  ;;
  g)
   COMPILER="gcc"
   COMPILER2="g++"
  ;;
  i)
   COMPILER="icc"
   COMPILER2="icc"
  ;;
  t)
   SIZE=""
  ;;
esac
done
if [ "`uname`" == "Darwin" ]; then
  PLATFORM="-D pp_OSX"
fi
export COMPILER
export SIZE
export PLATFORM
