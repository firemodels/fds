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
  else
    echo 1
    exit
  fi
}

gfortran_installed=`is_file_installed gfortran`
if [ $gfortran_installed -eq 0 ]; then
  echo unknown
  exit
fi

gfortran --version | head -1> gfortran_version 2>&1 
GFORTRANVERSION=`cat gfortran_version | awk '{print $4}' `
rm gfortran_version
echo "\"Gnu gfortran $GFORTRANVERSION\""
