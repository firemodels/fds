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
  if [ "$notfound" == "1" ] ; then
    echo 0
    exit
  fi
  echo 1
  exit
}

ifort_installed=`is_file_installed ifort`
if [ "$ifort_installed" == "0" ]; then
  echo Unknown
  exit
fi

ifort -v > ifort_version 2>&1 
IFORTVERSION=`cat ifort_version | awk '{print $3}' `
rm ifort_version
echo "\"Intel C/C++ $IFORTVERSION\""
