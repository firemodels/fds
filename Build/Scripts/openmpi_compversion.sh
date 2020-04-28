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

ifort_installed=`is_file_installed mpifort`
if [ $ifort_installed -eq 0 ]; then
  echo unknown
  exit
fi

mpifort -v > ifort_version 2>&1 
IFORTVERSION=`cat ifort_version | awk '{print $3}' `
rm ifort_version
echo "\"Intel ifort $IFORTVERSION\""
