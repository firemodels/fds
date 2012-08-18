#!/bin/csh -f

set target=$1
set fdsdir=$2
set host=$3
echo Building FDS in $fdsdir on $host using the target $target
ssh $host \( cd \~/$fdsdir \; rm -f *.o *.mod \; ./make_fds.sh \)
echo Completed
