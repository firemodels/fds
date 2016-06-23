#!/bin/csh -f
set dir=$1
set command=$2
set arg=$3
cd $dir
$command $arg
