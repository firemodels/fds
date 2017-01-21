#!/bin/bash
filein=$1
fileout=$filein.md5

MD5=md5sum
if [ "`uname`" == "Darwin" ] ; then
  MD5=md5
fi
notfound=`$MD5 -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ]; then
  exit
fi
if [ -e $filein ]; then
  $MD5 $filein > $fileout
fi
