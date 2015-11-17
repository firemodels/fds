#!/bin/bash
dir=`pwd`
target=${dir##*/}


echo Building $target
make VPATH="../../FDS_Source" -f ../makefile $target
