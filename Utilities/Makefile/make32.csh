#!/bin/csh -f
#
# set up environment for 32 bit Intel compilers
#
set target=$1
set DIR=/opt/intel/Compiler/11.0/074/bin/
source $DIR/iccvars.csh ia32
source $DIR/ifortvars.csh ia32
#
# build software
#
make $target
