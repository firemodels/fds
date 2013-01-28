#!/bin/bash
SIZE="-m64"
COMPILER="icc"
PLATFORM=""
while getopts '36ghio' OPTION
do
case $OPTION in
  h)
  echo "options:"
  echo "-3 - build a 32 bit version"
  echo "-6 - build a 64 bit version"
  echo "-g - use the gcc compiler"
  echo "-i - use the Intel icc compiler"
  echo "-o - built on an OSX platform"
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
  ;;
  i)
   COMPILER="icc"
  ;;
  o)
   PLATFORM="-D pp_OSX"
  ;;
esac
done
export COMPILER
export SIZE
export PLATFORM
