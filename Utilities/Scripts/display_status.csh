#!/bin/csh -f

set dir=$1
set file=$2.err
cd $dir
echo 
ls -l $file
tail -1 $file
