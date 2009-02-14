#!/bin/csh -f

cd $1

make VPATH="../../../FDS_Source" -f ../makefile intel_linux_32
