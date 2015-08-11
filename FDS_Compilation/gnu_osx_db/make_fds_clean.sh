#!/bin/bash
dir=`pwd`
target=${dir##*/}

echo Cleaning $target
make VPATH="../../FDS_Source" -f ../makefile clean
rm fds_$target
