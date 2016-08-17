#!/bin/bash

is_file_installed()
{
  program=$1
  notfound=`$program -help | tail -1 | grep "not found" | wc -l`
  if [ "$notfound" == "1" ] ; then
    echo "***error: $program not installed" 
    echo "          jpeg 2000 conversion aborted"
    exit
  fi
}

is_file_installed convert

shopt -s nullglob
for filename in *.jp2
do
  filebase=`basename $filename .jp2`
  echo converting $filebase.jp2 from jpeg 2000 to jpeg
  convert $filebase.jp2 $filebase.jpg
done
