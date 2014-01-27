#!/bin/csh -f
set dir=`pwd`
set target=$dir:t

echo Building $target
make VPATH="../../FDS_Source" -f ../makefile $target
