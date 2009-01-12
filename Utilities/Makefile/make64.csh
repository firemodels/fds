#!/bin/csh -f
#
# set up environment for 64 bit Intel compilers
#
set target=$1
set DIR=/opt/intel/Compiler/11.0/074/bin/
source $DIR/iccvars.csh intel64
source $DIR/ifortvars.csh intel64
#
# build software
#
make $target
