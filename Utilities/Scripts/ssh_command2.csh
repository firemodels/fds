#!/bin/csh -f
set host=$1
set dir=$2
set command=$3
set arg=$4
ssh -q $host \(cd $dir\; ./$command $arg\)
