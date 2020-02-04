#!/bin/bash

#---------------------------------------------
#                   is_file_installed
#---------------------------------------------

is_file_installed()
{
  local program=$1

  $program -help > prog_version 2>&1
  notfound=`cat prog_version | head -1 | grep "not found" | wc -l`
  rm prog_version
  if [ $notfound -eq 1 ] ; then
    echo 0
    exit
  fi
  echo 1
  exit
}

ifort_installed=`is_file_installed ifort`
if [ $ifort_installed -eq 0 ]; then
  echo "0"
  exit
fi

ifort -v >& version.out  
MAJORVERSION=`cat version.out | awk '{print $3}' | awk -F'.' '{print $1}'`
MINORVERSION=`cat version.out | awk '{print $3}' | awk -F'.' '{print $2}'`
rm version.out
if [[ "$MAJORVERSION" == "19" ]] && [[ "$MINORVERSION" != "0" ]]; then
  echo "20"
else
  echo "\"$MAJORVERSION\""
fi
