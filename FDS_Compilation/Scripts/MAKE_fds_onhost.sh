#!/bin/bash

target=$1
fdsdir=$2
host=$3
echo Building FDS in $fdsdir on $host using the target $target
ssh $host \( cd \~/$fdsdir \; rm -f *.o *.mod \; ./make_fds.sh \)
echo Completed
