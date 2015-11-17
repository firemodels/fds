#!/bin/bash
host=$1
dir=$2
command=$3
arg=$4
echo copying $dir/$command to $host
scp $dir/$command $host\:.
echo running command $command $arg on $host
ssh $host ./$command $arg
