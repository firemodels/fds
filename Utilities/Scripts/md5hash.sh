#!/bin/bash
filein=$1
fileout=$filein.md5

notfound=`md5sum -help 2>&1 | tail -1 | grep "not found" | wc -l`
if [ "$notfound" == "1" ]; then
  exit
fi
if [ -e $filein ]; then
  md5sum $filein > $fileout
fi
