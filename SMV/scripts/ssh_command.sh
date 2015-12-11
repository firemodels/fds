#!/bin/sh
host=$1
dir=$2
command=$3
arg=$4
#ssh -q $host bash -lc $dir/$command $arg
ssh -q $host $dir/$command $arg
