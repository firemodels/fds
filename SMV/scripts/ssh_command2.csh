#!/bin/csh -f
set host=$1
set dir=$2
set command=$3
set arg=$4
echo copying $dir/$command to $host
scp $dir/$command $host\:.
echo running command $command $arg on $host
ssh -q $host ./$command $arg
